"""AlphaFold3 docking wrapper.

Runs AF3 in no-MSA / single-sequence mode via a subprocess in a dedicated
conda env, then extracts ligand SDFs from the predicted mmCIFs using the
known input SMILES for bond-order recovery.

Public surface (used by cogligandbench.engine):
  - run_single(protein, ligand, output_dir, config=None, prefix=None, **kwargs)
  - run_dataset(config)
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional


def _smiles_from_sdf(sdf_path: str) -> str:
    """Read the first molecule from an SDF and return its canonical SMILES."""
    from rdkit import Chem

    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=True)
    except OSError as exc:
        raise ValueError(f"Cannot open SDF file {sdf_path}: {exc}") from exc
    mol = next(iter(suppl), None)
    if mol is None:
        raise ValueError(f"Cannot read molecule from {sdf_path}")
    return Chem.MolToSmiles(mol)


def _build_af3_input_json(system_id: str, pdb_path: str, sdf_path: str) -> Dict:
    """Construct the AF3 input-JSON dict for a single protein + ligand pair.

    Uses no-MSA mode: each protein chain's ``unpairedMsa`` field carries
    only the query sequence as a single-row A3M, ``pairedMsa`` is empty,
    and ``templates`` is the empty list.
    """
    from cogligandbench.utils.sequence import extract_protein_sequence

    sequences = extract_protein_sequence(pdb_path)
    if not sequences:
        raise ValueError(f"No protein sequences found in {pdb_path}")

    smiles = _smiles_from_sdf(sdf_path)

    entries: List[Dict] = []
    for i, seq in enumerate(sequences):
        entries.append({
            "protein": {
                "id": chr(ord("A") + i),
                "sequence": seq,
                "unpairedMsa": f">query\n{seq}\n",
                "pairedMsa": "",
                "templates": [],
            }
        })
    entries.append({"ligand": {"id": "L", "smiles": smiles}})

    return {
        "name": system_id,
        "modelSeeds": [1234],
        "sequences": entries,
        "dialect": "alphafold3",
        "version": 2,
    }


def _extract_ligand_from_cif(
    cif_path: str | Path,
    template_mol: "Chem.Mol",
) -> "Chem.Mol":
    """Pull the ligand (chain L) out of an AF3 mmCIF and recover bond orders.

    Strategy: parse the mmCIF with Biopython's MMCIF2Dict (a simple tokenizer
    that works with any minimal mmCIF loop), collect all atoms on chain L (the
    AF3 ligand chain) into a synthetic PDB block, parse that block with RDKit
    (which infers connectivity from distances), then call
    AssignBondOrdersFromTemplate with the known SMILES template to recover the
    correct bond orders.

    AF3 may label ligand atoms as either ATOM or HETATM depending on the
    upstream version; all atoms on chain L are collected regardless of the
    group_PDB label.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mmcif_dict = MMCIF2Dict(str(cif_path))

    # Extract per-atom arrays; chain id may live in auth_asym_id or label_asym_id
    chain_list = mmcif_dict.get(
        "_atom_site.auth_asym_id",
        mmcif_dict.get("_atom_site.label_asym_id", []),
    )
    resname_list = mmcif_dict.get("_atom_site.label_comp_id", [])
    atom_name_list = mmcif_dict.get("_atom_site.label_atom_id", [])
    elem_list = mmcif_dict.get("_atom_site.type_symbol", [])
    x_list = mmcif_dict.get("_atom_site.Cartn_x", [])
    y_list = mmcif_dict.get("_atom_site.Cartn_y", [])
    z_list = mmcif_dict.get("_atom_site.Cartn_z", [])

    n_atoms = len(x_list)
    if len(y_list) != n_atoms or len(z_list) != n_atoms or len(chain_list) != n_atoms:
        raise ValueError(
            f"Malformed mmCIF {cif_path}: _atom_site columns have inconsistent lengths "
            f"(x={len(x_list)}, y={len(y_list)}, z={len(z_list)}, chain={len(chain_list)})"
        )

    pdb_lines = []
    atom_idx = 1
    for i in range(n_atoms):
        if chain_list[i] != "L":
            continue
        resname = (resname_list[i] if i < len(resname_list) else "LIG")[:3].ljust(3)
        name = (atom_name_list[i] if i < len(atom_name_list) else "X")[:4].ljust(4)
        elem = (elem_list[i] if i < len(elem_list) else name.strip()[0]).strip()
        x = float(x_list[i])
        y = float(y_list[i])
        z = float(z_list[i])
        pdb_lines.append(
            f"HETATM{atom_idx:>5} {name} {resname} L{1:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          "
            f"{elem.rjust(2)}"
        )
        atom_idx += 1

    if not pdb_lines:
        raise ValueError(f"No ligand atoms (chain L) found in {cif_path}")

    pdb_block = "\n".join(pdb_lines) + "\nEND\n"
    raw_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
    if raw_mol is None:
        raise ValueError(f"Failed to parse synthesized PDB block for ligand from {cif_path}")

    try:
        Chem.SanitizeMol(raw_mol)
    except Exception:
        # Sanitization may fail before bond orders are recovered; that's OK.
        pass

    annotated = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
    return annotated


def _extract_ranked_ligand_sdfs(
    af3_system_dir: str | Path,
    smiles: str,
    out_dir: str | Path,
    num_poses: int,
) -> int:
    """Read ``ranking_scores.csv``, sort by score desc, write top-N rank{i}.sdf.

    Walks the per-system AF3 output directory, reads the ranking CSV, and
    writes the top ``num_poses`` ligand poses as ``rank{1..N}.sdf`` in
    ``out_dir``. Bond orders are recovered from the input ``smiles`` template.

    Returns the number of SDFs successfully written.
    """
    import pandas as pd
    from rdkit import Chem

    af3_system_dir = Path(af3_system_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scores_path = af3_system_dir / "ranking_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"AF3 ranking_scores.csv not found at {scores_path}")

    scores = pd.read_csv(scores_path)
    scores = scores.sort_values("ranking_score", ascending=False).reset_index(drop=True)

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"Could not parse template SMILES: {smiles}")

    written = 0
    for rank, row in scores.head(num_poses).iterrows():
        seed = int(row["seed"])
        sample = int(row["sample"])
        cif = af3_system_dir / f"seed-{seed}_sample-{sample}" / "model.cif"
        if not cif.exists():
            continue
        try:
            mol = _extract_ligand_from_cif(cif, template_mol=template)
        except Exception:
            continue
        Chem.MolToMolFile(mol, str(out_dir / f"rank{rank + 1}.sdf"))
        written += 1
    return written


def run_single(
    protein: str,
    ligand: str,
    output_dir: str,
    config: Optional[dict] = None,
    prefix: Optional[str] = None,
    **kwargs,
) -> str:
    """Dock a single protein+ligand pair with AlphaFold3.

    Implementation lands in a later task; this stub exists so the engine
    registration tests pass.
    """
    raise NotImplementedError("alphafold3 run_single not implemented yet")


def run_dataset(config: dict) -> None:
    """Run AlphaFold3 over an entire dataset directory.

    Implementation lands in a later task; this stub exists so the engine
    registration tests pass.
    """
    raise NotImplementedError("alphafold3 run_dataset not implemented yet")

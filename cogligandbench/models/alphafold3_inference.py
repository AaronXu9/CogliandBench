"""AlphaFold3 docking wrapper.

Runs AF3 in no-MSA / single-sequence mode via a subprocess in a dedicated
conda env, then extracts ligand SDFs from the predicted mmCIFs using the
known input SMILES for bond-order recovery.

Public surface (used by cogligandbench.engine):
  - run_single(protein, ligand, output_dir, config=None, prefix=None, **kwargs)
  - run_dataset(config)
"""

from __future__ import annotations

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


def _extract_ligand_from_cif(cif_path, template_mol):
    """Pull the ligand (chain L) out of an AF3 mmCIF and recover bond orders.

    Strategy: parse the mmCIF with Biopython's MMCIF2Dict (a simple tokenizer
    that works with any minimal mmCIF loop), collect HETATM records on chain
    "L" into a synthetic PDB block, parse that block with RDKit (which infers
    connectivity from distances), then call AssignBondOrdersFromTemplate with
    the known SMILES template to recover the correct bond orders.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mmcif_dict = MMCIF2Dict(str(cif_path))

    # Extract per-atom arrays; chain id may live in auth_asym_id or label_asym_id
    group_list = mmcif_dict.get("_atom_site.group_PDB", [])
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

    pdb_lines = []
    atom_idx = 1
    n_atoms = len(x_list)
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
        raise ValueError(f"No ligand atoms (chain L HETATM) found in {cif_path}")

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

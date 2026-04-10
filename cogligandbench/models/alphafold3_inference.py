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

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=True)
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

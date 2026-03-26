import glob
import os
import subprocess
import sys
import time
import traceback
import logging
from pathlib import Path

from omegaconf import OmegaConf
from rdkit import Chem
import rootutils

rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)
from cogligandbench.utils.log import get_custom_logger

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", str(rootutils.find_root(search_from=__file__, indicator=".project-root")))

AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HSD": "H", "HSE": "H", "HSP": "H", "HIE": "H", "HID": "H",
    "MSE": "M", "SEC": "U",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_protein_sequence(pdb_path: str) -> list:
    """Return a list of chain sequences from a PDB file."""
    from biopandas.pdb import PandasPdb
    ppdb = PandasPdb().read_pdb(pdb_path)
    atoms = ppdb.df["ATOM"]
    ca = atoms[atoms["atom_name"] == "CA"].drop_duplicates(
        subset=["chain_id", "residue_number", "insertion"]
    )
    sequences = []
    for chain_id, group in ca.groupby("chain_id", sort=False):
        seq = "".join(AA_3TO1.get(r, "X") for r in group["residue_name"])
        if seq:
            sequences.append(seq)
    return sequences


def _get_smiles_from_sdf(sdf_path: str) -> str:
    mol = next(Chem.SDMolSupplier(sdf_path, removeHs=True), None)
    if mol is None:
        raise ValueError(f"Cannot read molecule from {sdf_path}")
    return Chem.MolToSmiles(mol)


def _write_fasta(fasta_path: str, system_id: str, sequences: list, smiles: str):
    """Write a chai-compatible FASTA file."""
    with open(fasta_path, "w") as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">protein|{system_id}-chain-{i}\n{seq}\n")
        f.write(f">ligand|{system_id}-chain-{len(sequences) + 1}\n{smiles}\n")


def _run_chai_on_fasta(fasta_path: str, output_dir: str, config: dict):
    """Run Chai-1 inference on a FASTA file via subprocess (uses chai conda env python)."""
    python_exec = config.get("python_exec_path", "python3")
    cuda_device = config.get("cuda_device_index", 0)
    os.makedirs(output_dir, exist_ok=True)

    script = (
        "import torch\n"
        "from chai_lab.chai1 import run_inference\n"
        "from pathlib import Path\n"
        f"run_inference(\n"
        f"    fasta_file=Path('{fasta_path}'),\n"
        f"    output_dir=Path('{output_dir}'),\n"
        f"    num_trunk_recycles=3,\n"
        f"    num_diffn_timesteps=200,\n"
        f"    seed=42,\n"
        f"    device=torch.device('cuda:{cuda_device}'),\n"
        f"    use_esm_embeddings=True,\n"
        f")\n"
    )
    subprocess.run([python_exec, "-c", script], check=True)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_dataset(config: dict):
    """Run Chai-1 over a full dataset of FASTA input directories."""
    input_dir = config["input_dir"]
    output_dir = config["output_dir"]
    skip_existing = config.get("skip_existing", True)
    max_num_inputs = config.get("max_num_inputs", None)

    os.makedirs(output_dir, exist_ok=True)
    logger = get_custom_logger(
        "chai-1", config,
        f"chai-1_timing_{config.get('dataset', 'unknown')}_{config.get('repeat_index', 0)}.log",
    )

    num_processed = 0
    for item in sorted(os.listdir(input_dir)):
        item_path = os.path.join(input_dir, item)
        if not os.path.isdir(item_path):
            continue
        if max_num_inputs and num_processed >= max_num_inputs:
            logger.info(f"Reached max_num_inputs={max_num_inputs}. Stopping.")
            break

        out_item = os.path.join(output_dir, item)
        if skip_existing and os.path.exists(os.path.join(out_item, "pred.model_idx_0.pdb")):
            logger.info(f"Skipping {item} — output already exists.")
            continue

        fasta_files = glob.glob(os.path.join(item_path, "*.fasta"))
        if not fasta_files:
            logger.error(f"No FASTA found for {item}. Skipping.")
            continue

        os.makedirs(out_item, exist_ok=True)
        start = time.time()
        try:
            _run_chai_on_fasta(fasta_files[0], out_item, config)
            logger.info(f"{item},{time.time() - start:.2f}")
            num_processed += 1
        except Exception as e:
            logger.error(f"Failed {item}: {e}")
            with open(os.path.join(out_item, "error_log.txt"), "w") as f:
                traceback.print_exc(file=f)


def run_single(
    protein: str,
    ligand: str,
    output_dir: str,
    config: dict = None,
    prefix: str = None,
    **kwargs,
) -> str:
    """Run Chai-1 on a single protein (PDB) + ligand (SDF) pair.

    Generates a FASTA file and runs Chai-1 inference.
    Output: output_dir/pred.model_idx_*.pdb (Chai native format).
    """
    config = config or {}
    prefix = prefix or os.path.splitext(os.path.basename(protein))[0]
    os.makedirs(output_dir, exist_ok=True)

    # Build FASTA from PDB + SDF
    sequences = _extract_protein_sequence(protein)
    if not sequences:
        raise ValueError(f"No protein sequences found in {protein}")
    smiles = _get_smiles_from_sdf(ligand)

    fasta_path = os.path.join(output_dir, f"{prefix}.fasta")
    _write_fasta(fasta_path, prefix, sequences, smiles)

    try:
        _run_chai_on_fasta(fasta_path, output_dir, config)
        print(f"[INFO] Chai-1 docking complete. Results in {output_dir}")
    except Exception as e:
        print(f"[ERROR] Chai-1 docking failed: {e}")
        raise
    return output_dir


# ---------------------------------------------------------------------------
# Config loading + CLI entry point
# ---------------------------------------------------------------------------

def parse_config_file_with_omegaconf(config_path: str) -> dict:
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")
    os.environ.setdefault("PROJECT_ROOT", PROJECT_ROOT)
    cfg = OmegaConf.load(config_path)
    return OmegaConf.to_container(cfg, resolve=True)


if __name__ == "__main__":
    config = parse_config_file_with_omegaconf(sys.argv[1])
    run_dataset(config)

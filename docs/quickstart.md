# Quickstart — CogLigandBench

This guide covers getting started when you received the repository as a **directory copy** (rsync/scp), which includes the `data/` benchmark datasets. For a fresh GitHub clone, skip to [GitHub clone setup](#github-clone-setup).

---

## Directory copy setup

### 1. Install the Python package

From the repo root:

```bash
cd /path/to/CogLigandBench

# Create and activate a base conda env (Python 3.10 recommended)
conda create -n cogligandbench python=3.10 -y
conda activate cogligandbench

pip install -e .
```

This installs `cogligandbench` and its lightweight dependencies (rdkit, omegaconf, biopandas, pandas, etc.).

---

### 2. Configure environment variables

Copy the template and fill in paths for the methods you want to use:

```bash
cp .env.example .env
# edit .env with your paths
```

Then load them in your shell (add to `~/.bashrc` or source before running):

```bash
set -a && source .env && set +a
```

---

### 3. Method-by-method setup

#### Vina
Requires `obabel` and `vina` on `$PATH`.

```bash
# Ubuntu/Debian
sudo apt-get install openbabel
# vina binary: https://github.com/ccsb-scripps/AutoDock-Vina/releases
# place on PATH or add to ~/.bashrc
```

Test:
```bash
obabel --version && vina --version
```

---

#### GNINA
The binary is already bundled at `forks/GNINA/gnina` — no extra install needed. Requires a GPU.

---

#### UniDock2
The conda env is **not** bundled in the directory copy. Install it:

```bash
conda create -n unidock2 -c conda-forge unidock -y
```

---

#### DynamicBind
The conda env is bundled at `forks/DynamicBind/DynamicBind/`. However, conda envs contain hardcoded absolute paths and **may not work** if the directory was copied to a different path than the original.

**Check if it works:**
```bash
forks/DynamicBind/DynamicBind/bin/python3 -c "import torch; print(torch.__version__)"
```

**If it fails**, reinstall the env in place:
```bash
# Remove the broken env, recreate at the same prefix
rm -rf forks/DynamicBind/DynamicBind/
conda env create -f environments/dynamicbind_environment.yaml \
    --prefix forks/DynamicBind/DynamicBind/
```

Then set in `.env`:
```
DYNAMICBIND_PYTHON=/path/to/CogLigandBench/forks/DynamicBind/DynamicBind/bin/python3
DYNAMICBIND_DIR=/path/to/CogLigandBench/forks/DynamicBind
```

---

#### Chai-1
The conda env is bundled at `forks/chai-lab/chai-lab/`. Same portability caveat as DynamicBind.

**Check if it works:**
```bash
forks/chai-lab/chai-lab/bin/python -c "import chai_lab; print('OK')"
```

**If it fails**, reinstall:
```bash
rm -rf forks/chai-lab/chai-lab/
conda env create -f environments/chai_lab_environment.yaml \
    --prefix forks/chai-lab/chai-lab/
```

Chai-1 does not use an env var — pass the Python path directly to `dock_engine`:
```python
dock_engine('chai', protein=..., ligand=..., output_dir=...,
            python_exec_path='/path/to/CogLigandBench/forks/chai-lab/chai-lab/bin/python')
```

---

#### SurfDock
The conda env is **not** bundled — it lives at the system level. Install it:

```bash
conda env create -f environments/surfdock_environment.yaml -n SurfDock
```

SurfDock also requires a separate local installation with MSMS/APBS/pdb2pqr tools and an editable ESM install. See [SurfDock setup details](running_methods.md) for the full procedure.

Set in `.env`:
```
SURFDOCK_DIR=/path/to/your/SurfDock/installation
SURFDOCK_PRECOMPUTED_ARRAYS=/path/to/precomputed/precomputed_arrays
```

---

### 4. Verify the setup

Run the fast smoke tests (no docking, ~1s):

```bash
conda activate cogligandbench
pytest tests/test_engine_smoke.py -k "not slow" -v
```

All 25 should pass. Then test the methods you've set up:

```bash
# Test a specific method end-to-end (uses data/runsNposes/8gkf fixture)
pytest tests/test_engine_smoke.py -m slow -k "vina" -v
pytest tests/test_engine_smoke.py -m slow -k "gnina" -v
pytest tests/test_engine_smoke.py -m slow -k "dynamicbind" -v
pytest tests/test_engine_smoke.py -m slow -k "chai" -v
pytest tests/test_engine_smoke.py -m slow -k "unidock2" -v
pytest tests/test_engine_smoke.py -m slow -k "surfdock" -v
```

---

### 5. Run docking

```python
from cogligandbench import dock_engine

# Single molecule
result = dock_engine(
    'vina',
    protein='my_receptor.pdb',
    ligand='my_ligand.sdf',
    output_dir='./results',
)

# Run over the full runsNposes benchmark
dock_engine('vina', dataset='runsNposes')
dock_engine('gnina', dataset='runsNposes', repeat_index=1)
```

Poses are written as ranked SDF files in `result/`. Compute RMSD against the crystal ligand:

```python
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from pathlib import Path

ref  = Chem.MolFromMolFile('crystal_ligand.sdf', removeHs=True)
pose = Chem.MolFromMolFile(str(next(Path(result).glob('*rank1*.sdf'))), removeHs=True)
rmsd = CalcRMS(ref, pose)
print(f'RMSD: {rmsd:.2f} Å')
```

---

## GitHub clone setup

If you cloned from GitHub (no data included), the workflow is identical except you also need to download the benchmark data.

```bash
git clone https://github.com/AaronXu9/CogliandBench --recursive
cd CogliandBench
pip install -e .
```

The `data/` directory is not tracked in git. Contact the lab for access to the benchmark datasets, or download the publicly available ones:

| Dataset | Source |
|---------|--------|
| runsNposes | [runs-n-poses GitHub](https://github.com/plinder-org/runs-n-poses) |
| PoseBusters | [posebusters GitHub](https://github.com/maabuu/posebusters) |
| Astex Diverse | PDB (21 complexes) |
| DockGen | [DockGen GitHub](https://github.com/HannesStark/DockGen) |

After placing data under `data/<dataset_name>/`, follow the method setup in [Step 3](#3-method-by-method-setup) above.

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'cogligandbench'`**
→ Run `pip install -e .` from the repo root with your conda env active.

**`KeyError: 'SURFDOCK_DIR'` or similar**
→ Your `.env` is missing that variable, or you forgot to `source .env`. Check `.env.example` for the full list.

**Conda env Python import errors after copy**
→ The env was built with absolute paths to the original machine. Reinstall it at the same prefix as shown in the method setup above.

**GNINA runs on CPU only / very slow**
→ GNINA requires CUDA. Check `nvidia-smi` and ensure the GNINA binary was compiled with GPU support.

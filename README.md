# CogLigandBench

![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

A benchmarking framework for protein-ligand docking methods evaluated against experimental crystal structures.

---

## Why Crystal Structures

Most docking benchmarks pair a ligand with an AlphaFold-predicted or homology-modelled receptor. CogLigandBench uses **co-crystal structures from the PDB** instead: the receptor conformation is the one that was experimentally observed in the presence of the ligand. This isolates docking error from structure-prediction error and mirrors the real structure-based drug design (SBDD) workflow, where a solved crystal structure of the target is the starting point.

**Benchmark datasets:**

| Dataset | Systems | Source |
|---------|---------|--------|
| runsNposes | ~1,280 | [runs-n-poses](https://github.com/plinder-org/runs-n-poses) |
| Astex Diverse | 85 | PDB (curated fragment set) |
| PoseBusters | 428 | [posebusters](https://github.com/maabuu/posebusters) |
| DockGen | 189 | [DockGen](https://github.com/HannesStark/DockGen) |

---

## Supported Methods & Preliminary Results

Results below are **top-1 pose, runsNposes benchmark** (crystal structures, no structure prediction).

| Method | Type | % ≤ 2 Å RMSD | Median RMSD |
|--------|------|-------------|-------------|
| SurfDock | Deep learning | 59.5% | 1.56 Å |
| GNINA | CNN scoring | 33.9% | 2.56 Å |
| ICM-RTCNN | Physics + ML | 27.5% | 3.42 Å |
| ICM | Physics-based | 25.7% | 3.69 Å |
| Vina | Physics-based | 7.5% | 7.43 Å |
| UniDock2 | Physics-based | 2.0%† | 24.35 Å† |

†UniDock2 result is tentative — run artifact under investigation.

**Conda environments by method:**

| Method | Conda env | Notes |
|--------|-----------|-------|
| `vina` | system | Needs `obabel`, `vina` on `$PATH` |
| `gnina` | system | Binary at `forks/GNINA/gnina` |
| `chai` | `forks/chai-lab/chai-lab/` | Pass `python_exec_path=` to `dock_engine` |
| `dynamicbind` | `forks/DynamicBind/DynamicBind/` | Path set in YAML config |
| `unidock2` | `unidock2` | `conda create -n unidock2 -c conda-forge unidock` |
| `surfdock` | `SurfDock` | Requires MSMS/APBS/pdb2pqr tools |
| `icm` | system | Requires commercial ICM license |

---

## Installation & Quickstart

### 1. Install the package

```bash
cd /path/to/CogLigandBench
conda create -n cogligandbench python=3.10 -y
conda activate cogligandbench
pip install -e .
```

### 2. Configure environment variables

```bash
cp .env.example .env
# edit .env with paths for the methods you want to use
set -a && source .env && set +a
```

### 3. Run docking

```python
from cogligandbench import dock_engine

# Single molecule
result = dock_engine(
    'vina',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./results',
)

# Full runsNposes benchmark
dock_engine('gnina', dataset='runsNposes')
dock_engine('gnina', dataset='runsNposes', repeat_index=1)
```

Poses are written as ranked SDF files under `results/`. See [docs/quickstart.md](docs/quickstart.md) for per-method setup and [docs/api.md](docs/api.md) for the full `dock_engine` API.

### 4. Compute RMSD

```python
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from pathlib import Path

ref  = Chem.MolFromMolFile('crystal_ligand.sdf', removeHs=True)
pose = Chem.MolFromMolFile(str(next(Path(result).glob('*rank1*.sdf'))), removeHs=True)
rmsd = CalcRMS(ref, pose)
print(f'RMSD: {rmsd:.2f} Angstrom')
```

---

## Project Structure

```
cogligandbench/         # Python package -- dock_engine, per-method wrappers
  engine.py             # dock_engine unified entry point
  models/               # per-method inference modules (run_dataset, run_single)
  data/                 # input preparation and output extraction scripts
  analysis/             # RMSD, complex alignment, scoring
  utils/                # logging and data utilities

cogligand_config/       # Hydra/OmegaConf YAML configs
  model/                # per-method inference configs

data/                   # benchmark datasets (crystal structures, not in git)
  runsNposes/
  astex_diverse_set/
  posebusters_benchmark_set/
  dockgen_set/

forks/                  # third-party codebases (submodules / local installs)
  GNINA/                # GNINA binary
  Vina/                 # AutoDock Vina + ADFR suite
  DynamicBind/          # DynamicBind diffusion docking
  chai-lab/             # Chai-1
  UniDock2/             # UniDock2
  SurfDock/             # SurfDock (code + weights; tools installed separately)
  ICM/                  # ICM docking scripts

docs/                   # documentation
tests/                  # pytest suite (fast smoke tests + slow end-to-end)
```

---

## Citation & Acknowledgements

If you use CogLigandBench, please cite:

```bibtex
@software{cogligandbench2024,
  author  = {Xu, Aaron and others},
  title   = {CogLigandBench: Crystal-Structure Protein-Ligand Docking Benchmark},
  year    = {2024},
  url     = {https://github.com/AaronXu9/CogliandBench},
}
```

CogLigandBench was built on top of [PoseBench](https://github.com/BioinfoMachineLearning/PoseBench) (Morehead et al., 2024), which targets predicted/AlphaFold structures.

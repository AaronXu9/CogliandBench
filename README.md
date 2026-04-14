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
| runsNposes | ~1,280 | [runs-n-poses](https://github.com/plinder-org/runs-n-poses) — primary benchmark |
| Plinder test set | ~1,038 | [plinder](https://www.plinder.sh/) — secondary validation |

---

## Supported Methods & Benchmark Results

Results below are RMSD < 2 Å success rates on the **runsNposes** benchmark, stratified by SuCOS similarity to pre-2021 training data. Source: paper Table 1 (Xu et al., *A Multi-Paradigm Benchmark of Molecular Docking*, see `docs/A_Multi-Paradigm_Benchmark_main_v8.docx`).

| Category | Method | Success 0–20% (%) | Success 80–100% (%) | Gap Δ (pp) |
|----------|--------|-------------------|---------------------|------------|
| Reference  | Best Possible | 89 | 100 | +11 |
| Co-folding | AlphaFold3    | 38 |  96 | +58 |
| Co-folding | Boltz-1       | 22 |  94 | +72 |
| Co-folding | Chai-1        | 25 |  94 | +69 |
| Co-folding | Protenix      | 27 |  92 | +65 |
| Hybrid     | GNINA         | 62 |  89 | +27 |
| Hybrid     | ICM-RTCNN     | 67 |  79 | +12 |
| Physics    | ICM           | 53 |  71 | +18 |
| Physics    | Vina          | 18 |  31 | +13 |
| Physics    | Uni-Dock      | 12 |  27 | +15 |

The **generalization gap (Δ)** is the success-rate difference between the highest- and lowest-similarity bins. Co-folding models show 58–72 pp gaps (memorization-sensitive); hybrid and physics methods are more OOD-robust. Boltz-2 is implemented in this repo but excluded from the table because its training set overlaps the test set.

### Pose Comparison: Co-folding vs Hybrid in Two Similarity Regimes

Predicted poses (magenta sticks) overlaid on the crystal ligand (green sticks) and receptor (grey cartoon) for one **high-similarity** target (`7qhl__1__1.A__1.E`, SuCOS = 85.2) and one **low-similarity** target (`8cfb__1__1.A__1.L`, SuCOS = 3.0). AF3's predicted protein is aligned to the crystal receptor before showing only its ligand.

| | High similarity (SuCOS 85.2) | Low similarity (SuCOS 3.0) |
|---|---|---|
| **AF3** (co-folding) | ![AF3 high sim](docs/figures/af3_high_sim.png) <br> RMSD = 0.19 Å (success) | ![AF3 low sim](docs/figures/af3_low_sim.png) <br> RMSD = 37.8 Å (hallucinated, far from pocket) |
| **GNINA** (hybrid) | ![GNINA high sim](docs/figures/gnina_high_sim.png) <br> RMSD = 1.42 Å (success) | ![GNINA low sim](docs/figures/gnina_low_sim.png) <br> RMSD = 0.44 Å (success) |

At high similarity both paradigms recover the native pose. At low similarity AF3 places the ligand far outside the binding site (RMSD ≈ 38 Å) while GNINA's physics-based sampling still locks into the correct pose (RMSD < 0.5 Å). This is the contrast quantified by the generalization gap in Table 1.

The figures were generated with `scripts/generate_pose_comparison_figure.py`; see `docs/figures/selection_metadata.json` for the exact selection criteria.

---

**Conda environments by method:**

Each supported method has an idempotent install script under `scripts/install_{method}_env.sh`. Running it creates the env at `/mnt/katritch_lab2/aoxu/envs/{method}` and symlinks `envs/{method}` from the project root.

| Method | Install script | Env path (after install) |
|--------|----------------|--------------------------|
| `alphafold3` | `scripts/install_alphafold3_env.sh` | `envs/alphafold3/` |
| `protenix`   | `scripts/install_protenix_env.sh`   | `envs/protenix/` |
| `boltz1`/`boltz2` | `scripts/install_boltz_env.sh` | `envs/boltz/` |
| `chai`       | `scripts/install_chai_env.sh`       | `envs/chai/` |
| `vina`       | `scripts/install_vina_env.sh`       | `envs/vina/` |
| `gnina`      | `scripts/install_gnina_env.sh`      | `envs/gnina/` (wraps `forks/GNINA/`) |
| `dynamicbind`| `scripts/install_dynamicbind_env.sh`| `envs/dynamicbind/` (wraps `forks/DynamicBind/`) |
| `unidock2`   | `scripts/install_unidock2_env.sh`   | `envs/unidock2/` |
| `surfdock`   | `scripts/install_surfdock_env.sh`   | `envs/surfdock/` |
| `icm`        | — | physics-based; uses commercial ICM binary on `$PATH` |

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
  runsNposes/           # primary benchmark
  plinder_set/          # secondary validation

envs/                   # symlinks to per-method conda envs (created by install scripts)
  alphafold3/, protenix/, boltz/, chai/, vina/, gnina/, ...

forks/                  # third-party codebases (submodules / local installs)
  alphafold3/           # AlphaFold3 source + weights
  GNINA/                # GNINA binary
  Vina/                 # AutoDock Vina + ADFR suite
  DynamicBind/          # DynamicBind diffusion docking
  chai-lab/             # Chai-1 source
  UniDock2/             # UniDock2
  SurfDock/             # SurfDock (code + weights; tools installed separately)
  ICM/                  # ICM docking scripts (commercial binary required)
  boltz/                # Boltz support files (input prep, extraction)

docs/                   # documentation
tests/                  # pytest suite (fast smoke tests + slow end-to-end)
```

---

## Citation

If you use CogLigandBench, please cite:

```bibtex
@article{xu2026multiparadigm,
  author  = {Xu, Ao and Lam, Jordy Homing and Nakano, Aiichiro and Katritch, Vsevolod},
  title   = {A Multi-Paradigm Benchmark of Molecular Docking: From Physics to Co-folding and Hybrid Models},
  year    = {2026},
  note    = {Manuscript; data and code at \url{https://zenodo.org/records/16754298} and \url{https://github.com/AaronXu9/CogligandBench}},
}
```

The codebase was originally seeded from [PoseBench](https://github.com/BioinfoMachineLearning/PoseBench) (Morehead et al., 2024) and has since been extensively rewritten for crystal-structure benchmarks.

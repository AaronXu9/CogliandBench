# API Reference — CogLigandBench

## `dock_engine`

```python
from cogligandbench import dock_engine

result = dock_engine(
    method,
    *,
    # dataset mode
    dataset=None,
    repeat_index=0,
    # single-molecule mode
    protein=None,
    ligand=None,
    output_dir=None,
    prefix=None,
    # per-method kwargs (merged into YAML config)
    **kwargs,
)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | `str` | required | One of `'vina'`, `'gnina'`, `'chai'`, `'dynamicbind'`, `'unidock2'`, `'surfdock'` |
| `dataset` | `str | None` | `None` | Dataset name for batch mode (e.g. `'runsNposes'`). Mutually exclusive with `protein`/`ligand`. |
| `repeat_index` | `int` | `0` | Run index; used in output directory naming for batch mode. |
| `protein` | `str | Path | None` | `None` | Path to receptor PDB for single-molecule mode. |
| `ligand` | `str | Path | None` | `None` | Path to ligand SDF for single-molecule mode. |
| `output_dir` | `str | Path | None` | `None` | Output root for single-molecule mode. Results land in `output_dir/{method}/{prefix}/`. |
| `prefix` | `str | None` | stem of protein filename | Identifier used for output filenames in single-molecule mode. |
| `**kwargs` | any | — | Method-specific overrides merged into the YAML config (see per-method tables below). |

### Return value

- **Single-molecule mode:** `str` — path to `output_dir/{method}/{prefix}/` containing ranked pose files.
- **Dataset mode:** `None` — results written directly to the path in the YAML config.

---

## Single-molecule vs dataset mode

```python
# Single-molecule mode — provide protein, ligand, output_dir
result_dir = dock_engine(
    'vina',
    protein='receptor.pdb',
    ligand='ligand.sdf',
    output_dir='./results',
    prefix='my_system',
)
# Poses written to ./results/vina/my_system/

# Dataset mode — provide dataset name
dock_engine('vina', dataset='runsNposes')
dock_engine('vina', dataset='runsNposes', repeat_index=1)
# Poses written to path defined in cogligand_config/model/vina_inference.yaml
```

---

## Per-method kwargs

All methods accept `**kwargs` that override keys in their YAML config.

### Vina

| Kwarg | Type | Default | Description |
|-------|------|---------|-------------|
| `top_n` | `int` | `10` | Number of poses to write |
| `exhaustiveness` | `int` | `32` | Vina search exhaustiveness |
| `num_modes` | `int` | `20` | Max number of binding modes |
| `box_size` | `list[float]` | auto | `[x, y, z]` in Å; auto-computed from ligand if omitted |

### GNINA

| Kwarg | Type | Default | Description |
|-------|------|---------|-------------|
| `top_n` | `int` | `10` | Number of poses to write |
| `exhaustiveness` | `int` | `32` | Search exhaustiveness |
| `num_modes` | `int` | `20` | Max binding modes |
| `cnn_scoring` | `str` | `'rescore'` | CNNscore mode: `'rescore'`, `'refinement'`, `'metrorescore'` |

### Chai-1

| Kwarg | Type | Default | Description |
|-------|------|---------|-------------|
| `python_exec_path` | `str` | `'python3'` | **Required in practice.** Path to the chai conda env Python. |
| `num_trunk_recycles` | `int` | `3` | Number of trunk recycling iterations |
| `num_diffn_timesteps` | `int` | `200` | Diffusion timesteps |
| `seed` | `int` | `42` | Random seed |

### DynamicBind

| Kwarg | Type | Default | Description |
|-------|------|---------|-------------|
| `samples_per_complex` | `int` | `40` | Total samples to generate |
| `savings_per_complex` | `int` | `40` | Top-N samples to save |
| `python_exec_path` | `str` | from YAML | DynamicBind conda env Python (set in config) |

### UniDock2

| Kwarg | Type | Default | Description |
|-------|------|---------|-------------|
| `num_poses` | `int` | `10` | Number of output poses |
| `box_size` | `list[float]` | `[15, 15, 15]` | Docking box in Å |
| `search_mode` | `str` | `'detail'` | `'fast'`, `'balance'`, or `'detail'` |
| `energy_range` | `float` | `15.0` | Energy window in kcal/mol |

### SurfDock

| Kwarg | Type | Default | Description |
|-------|------|---------|-------------|
| `num_poses` | `int` | `40` | Number of poses to keep |
| `samples_per_complex` | `int` | `40` | Diffusion samples to generate |
| `batch_size` | `int` | `5` | Inference batch size |
| `surfdock_dir` | `str` | `$SURFDOCK_DIR` | Path to local SurfDock installation |

---

## Output file naming

| Method | Pattern | Location |
|--------|---------|----------|
| `vina` | `{prefix}_pose{N}_score{S:.2f}.sdf` | `output_dir/vina/{prefix}/` |
| `gnina` | `{prefix}_pose{N}_score{S:.2f}.sdf` | `output_dir/gnina/{prefix}/` |
| `chai` | `pred.model_idx_{0-4}.pdb` | `output_dir/chai/{prefix}/` |
| `dynamicbind` | `rank{N}_ligand_lddt{L}_affinity{A}.sdf` | `output_dir/dynamicbind/{prefix}/` |
| `unidock2` | `rank{N}.sdf` | `output_dir/unidock2/{prefix}/` |
| `surfdock` | `rank{N}.sdf` | `output_dir/surfdock/{prefix}/` |

**Ranking:** Lower rank number = better pose, except GNINA where `score` is CNNscore (higher = better).

---

## Config override pattern

Each `**kwarg` passed to `dock_engine` is merged into the method's YAML config before inference. Config files live in `cogligand_config/model/<method>_inference.yaml`.

```python
# This call:
dock_engine('vina', protein='r.pdb', ligand='l.sdf', output_dir='./out', top_n=5)

# Is equivalent to loading cogligand_config/model/vina_inference.yaml
# then overriding: cfg.top_n = 5
```

To inspect the full list of config keys for a method:

```python
from omegaconf import OmegaConf
cfg = OmegaConf.load('cogligand_config/model/vina_inference.yaml')
print(OmegaConf.to_yaml(cfg))
```

---

## RMSD calculation

```python
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from pathlib import Path


def top1_rmsd(crystal_sdf: str, result_dir: str, method: str = 'vina') -> float:
    """Compute symmetry-corrected RMSD for the top-ranked pose."""
    ref = Chem.MolFromMolFile(crystal_sdf, removeHs=True)

    if method == 'chai':
        mol = Chem.MolFromPDBFile(
            str(next(Path(result_dir).glob('pred.model_idx_0.pdb'))),
            removeHs=True,
        )
    elif method in ('unidock2', 'surfdock', 'dynamicbind'):
        mol = Chem.MolFromMolFile(
            str(next(Path(result_dir).glob('rank1*.sdf'))),
            removeHs=True,
        )
    else:  # vina, gnina
        mol = Chem.MolFromMolFile(
            str(next(Path(result_dir).glob('*pose1*.sdf'))),
            removeHs=True,
        )

    return CalcRMS(ref, mol)  # symmetry-corrected
```

**Note for Chai-1:** The `pred.model_idx_*.pdb` files contain the full protein-ligand complex. Extract the ligand residue before computing RMSD, or use `cogligandbench/data/chai_output_extraction.py` for batch extraction.

# RMSD Analysis — CogLigandBench

This document describes how to compute RMSD between docked poses and crystal ligands for all supported methods. All analysis uses RDKit's symmetry-corrected RMSD (`CalcRMS`), which accounts for equivalent atom orderings.

---

## Output locations and file patterns

After running `dock_engine`, poses land in `output_dir/{method}/{prefix}/`. The table below lists the rank-1 file pattern for each method.

| Method | Rank-1 file | Notes |
|--------|------------|-------|
| `vina` | `*_pose1_score*.sdf` | Scored by Vina energy (most negative = best) |
| `gnina` | `*_pose1_score*.sdf` | Scored by CNNscore (highest = best) |
| `chai` | `pred.model_idx_0.pdb` | Full complex PDB; ligand must be extracted |
| `dynamicbind` | `rank1_ligand_lddt*.sdf` | Scored by DynamicBind lDDT |
| `unidock2` | `rank1.sdf` | Scored by Vina energy (most negative = best) |
| `surfdock` | `rank1.sdf` | Scored by confidence (highest = best) |

---

## Single-system RMSD

```python
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
from pathlib import Path


def compute_rmsd(crystal_sdf: str, pose_sdf: str) -> float:
    ref  = Chem.MolFromMolFile(crystal_sdf, removeHs=True)
    pose = Chem.MolFromMolFile(pose_sdf, removeHs=True)
    if ref is None or pose is None:
        raise ValueError(f"Could not load molecule from {crystal_sdf} or {pose_sdf}")
    return CalcRMS(ref, pose)  # symmetry-corrected


# Example
rmsd = compute_rmsd(
    'data/runsNposes/8gkf__1__1.A__1.J/8gkf__1__1.A__1.J_ligand.sdf',
    'results/vina/8gkf_test/8gkf_test_pose1_score-7.23.sdf',
)
print(f'RMSD: {rmsd:.2f} Angstrom')
```

---

## Batch RMSD over runsNposes

```python
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS

DATA_DIR   = Path('data/runsNposes')   # crystal structures
RESULTS    = Path('results/vina')      # dock_engine output root
RANK1_GLOB = '*_pose1_score*.sdf'      # adjust per method (see table above)

records = []
for system_dir in sorted(DATA_DIR.iterdir()):
    system = system_dir.name
    crystal = system_dir / f'{system}_ligand.sdf'
    pose_dir = RESULTS / system
    poses = sorted(pose_dir.glob(RANK1_GLOB)) if pose_dir.exists() else []

    if not crystal.exists() or not poses:
        records.append({'system': system, 'rmsd': None})
        continue

    ref  = Chem.MolFromMolFile(str(crystal), removeHs=True)
    pose = Chem.MolFromMolFile(str(poses[0]), removeHs=True)
    try:
        rmsd = CalcRMS(ref, pose)
    except Exception:
        rmsd = None
    records.append({'system': system, 'rmsd': rmsd})

df = pd.DataFrame(records)
valid = df['rmsd'].dropna()
print(f"Systems:         {len(df)}")
print(f"With results:    {valid.count()}")
print(f"Median RMSD:     {valid.median():.2f} Angstrom")
print(f"% <= 2 Ang (top-1): {(valid <= 2.0).mean() * 100:.1f}%")
print(f"% <= 5 Ang (top-1): {(valid <= 5.0).mean() * 100:.1f}%")
```

To run this for a different method, change `RESULTS` and `RANK1_GLOB`:

| Method | `RESULTS` | `RANK1_GLOB` |
|--------|-----------|-------------|
| gnina | `Path('results/gnina')` | `'*_pose1_score*.sdf'` |
| dynamicbind | `Path('results/dynamicbind')` | `'rank1_ligand_lddt*.sdf'` |
| unidock2 | `Path('results/unidock2')` | `'rank1.sdf'` |
| surfdock | `Path('results/surfdock')` | `'rank1.sdf'` |

---

## Chai-1: extract ligand before RMSD

Chai-1 outputs full protein-ligand complex PDBs. Extract the ligand first:

```python
from cogligandbench.data.chai_output_extraction import extract_ligand_from_complex

# Extract ligand from complex PDB
ligand_sdf = extract_ligand_from_complex(
    complex_pdb='results/chai/8gkf_test/pred.model_idx_0.pdb',
    output_sdf='results/chai/8gkf_test/ligand_model_0.sdf',
)

# Then compute RMSD normally
rmsd = compute_rmsd(
    'data/runsNposes/8gkf__1__1.A__1.J/8gkf__1__1.A__1.J_ligand.sdf',
    ligand_sdf,
)
```

For batch extraction over a full dataset run, use `cogligandbench/data/chai_output_extraction.py`.

---

## Preliminary benchmark results (runsNposes, top-1)

| Method | Systems | % <= 2 Ang | % <= 5 Ang | Median RMSD |
|--------|---------|-----------|-----------|-------------|
| SurfDock | ~1,280 | 59.5% | 79.2% | 1.56 Ang |
| GNINA | ~1,280 | 33.9% | 58.1% | 2.56 Ang |
| ICM-RTCNN | ~1,280 | 27.5% | 54.3% | 3.42 Ang |
| ICM | ~1,280 | 25.7% | 52.8% | 3.69 Ang |
| Vina | ~1,280 | 7.5% | 28.4% | 7.43 Ang |
| UniDock2 | ~1,280 | 2.0%+ | — | 24.35 Ang+ |

+UniDock2 tentative — run artifact under investigation.

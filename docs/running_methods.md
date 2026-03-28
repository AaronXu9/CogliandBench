# Running Docking Methods (cogligandbench)

All methods operate on **crystal structures** via the `cogligandbench` package. Scripts are run from the repo root. The `PROJECT_ROOT` environment variable must be set (or resolved via `rootutils`/`.project-root`).

---

## Common Pattern

Most methods follow the same invocation pattern:

```bash
export PROJECT_ROOT=/path/to/CogLigandBench

# Hydra-based methods (chai, dynamicbind):
python3 cogligandbench/models/<method>_inference.py [key=value overrides]

# Config-file-based methods (vina, gnina, ICM):
python3 <script_path> cogligand_config/model/<method>_inference.yaml
```

Config files are in `cogligand_config/model/`. Key fields common to all configs:
- `repeat_index` — run index (used in output dir naming)
- `benchmark` / `dataset` — dataset name
- `skip_existing` — resume incomplete runs
- `logging.log_dir` — where timing logs are written

---

## Vina

**Script:** `cogligandbench/models/vina_inference.py`
**Config:** `cogligand_config/model/vina_inference.yaml`

```bash
python3 cogligandbench/models/vina_inference.py cogligand_config/model/vina_inference.yaml
```

The config is loaded via OmegaConf (not Hydra). Key config fields:
```yaml
benchmark: runsNposes
repeat_index: 0
method: vina
inputs_csv: ${oc.env:PROJECT_ROOT}/forks/Vina/inference/vina_runsNposes_benchmark_inputs.csv
output_dir: ${oc.env:PROJECT_ROOT}/forks/Vina/inference/${benchmark}_${repeat_index}
log_dir: ${oc.env:PROJECT_ROOT}/forks/Vina/inference/logs/Vina
top_n: 10
```

The script reads `inputs_csv` (columns: `complex_name`, `protein_path`, `ligand_path`), converts protein/ligand to PDBQT via OpenBabel, computes a docking box from the reference ligand (with MDAnalysis-based refinement if available), runs `vina`, and parses top-N poses to SDF via OpenBabel/pybel.

To switch to GNINA scoring instead of Vina docking, set `method: gnina` in the config — this uses `forks/GNINA/gnina` with `--autobox_ligand`.

---

## ICM

ICM has two separate scripts depending on the use case:

### General datasets (`dock_refactor.py`)
**Script:** `forks/ICM/icm/dock_refactor.py`
**Config:** `cogligand_config/model/icm_inference.yaml`

```bash
python3 forks/ICM/icm/dock_refactor.py cogligand_config/model/icm_inference.yaml
```

### runsNposes dataset (`dock_runsNposes.py`)
**Script:** `forks/ICM/icm/dock_runsNposes.py`
**Config:** `cogligand_config/model/icm_inference.yaml`

```bash
python3 forks/ICM/icm/dock_runsNposes.py cogligand_config/model/icm_inference.yaml
```

Key config fields:
```yaml
dataset: runsNposes
task: dock          # or: identify_pocket, rank_results
data_dir: /path/to/CogLigandBench/data/runsNposes/
icm_executable: /home/aoxu/icm-3.9-4/icm64
icm_dockscan_path: /home/aoxu/icm-3.9-4/_dockScan
icm_docking_dir: /path/to/CogLigandBench/forks/ICM
icm_map_dir: /path/to/CogLigandBench/forks/ICM/ICM_manual_docking_maps/runsNposes
icb_out_dir: /path/to/CogLigandBench/forks/ICM/inference/runsNposes
docking_maps: manual
docking_params:
  num_conf: 10
  thorough: 10.0
  effort: 3
skip_existing: true
```

ICM tasks run sequentially: `identify_pocket` → `dock` → `rank_results`. Protein names longer than 24 characters are automatically shortened (a JSON mapping is saved to `icm_map_dir/protein_name_map.json`). ICM templates are in `forks/ICM/icm_docking_scripts_template/`.

Platform paths are auto-detected (Darwin vs Linux) via `platform_configs` in the config.

---

## Chai-1

**Script:** `cogligandbench/models/chai_inference.py`
**Config:** `cogligand_config/model/chai_inference.yaml`

```bash
python3 cogligandbench/models/chai_inference.py
# or with overrides:
python3 cogligandbench/models/chai_inference.py dataset=runsNposes repeat_index=0
```

Key config fields:
```yaml
dataset: runsNposes
input_dir: ${oc.env:PROJECT_ROOT}/forks/chai-lab/prediction_inputs/${dataset}
output_dir: ${oc.env:PROJECT_ROOT}/forks/chai-lab/prediction_outputs/${dataset}_${repeat_index}
cuda_device_index: 0
repeat_index: 0
skip_existing: true
pocket_only_baseline: false
```

Inputs are per-complex directories under `input_dir/`, each containing a `.fasta` file. Chai-1 is invoked via `chai_lab.chai1.run_inference` with `num_trunk_recycles=3`, `num_diffn_timesteps=200`, `seed=42`. Outputs land in `output_dir/<complex_name>/pred.model_idx_0.pdb`.

---

## DynamicBind

**Script:** `cogligandbench/models/dynamicbind_inference.py`
**Config:** `cogligand_config/model/dynamicbind_inference.yaml`

```bash
python3 cogligandbench/models/dynamicbind_inference.py dataset=runsNposes repeat_index=0
```

Reads protein files from `input_data_dir` and ligand CSV files from `input_ligand_csv_dir`. Pairs proteins and ligands by matching filename stems. Runs DynamicBind via subprocess from `forks/DynamicBind/`.

---

## SurfDock

**Script:** `cogligandbench/models/surfdock_inference.py`
**Config:** `cogligand_config/model/surfdock_inference.yaml`

SurfDock requires a local installation with MSMS/APBS/pdb2pqr binaries. Set the required env vars before running:

```bash
export SURFDOCK_DIR=/path/to/your/SurfDock/installation
export SURFDOCK_PRECOMPUTED_ARRAYS=/path/to/precomputed/precomputed_arrays
```

For dataset mode, pre-computed inputs are required:

```bash
export SURFDOCK_INPUTS_CSV=/path/to/surfdock_runsNposes_benchmark_inputs.csv
export SURFDOCK_ESM_EMBEDDINGS=/path/to/esm2_3billion_pdbbind_embeddings.pt
export SURFDOCK_SURFACE_DIR=/path/to/8A_surface
```

---

## Notes

- **Input CSV format** (Vina/GNINA): columns `complex_name`, `protein_path`, `ligand_path`
- **Input directory format** (Chai-1): `<input_dir>/<complex_name>/<complex_name>.fasta`
- **Input directory format** (ICM): `<data_dir>/<complex_name>/<complex_name>.pdb` + `*.sdf`
- **Output poses** are written as ranked SDF files: `<complex_name>_pose<rank>_score<score>.sdf`
- **Timing logs** record `complex_name,elapsed_seconds` per entry

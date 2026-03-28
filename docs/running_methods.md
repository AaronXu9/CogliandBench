# Running Docking Methods (cogligandbench)

All methods operate on **crystal structures** via the `cogligandbench` package. Scripts are run from the repo root. The `PROJECT_ROOT` environment variable must be set (or resolved via `rootutils`/`.project-root`).

---

## Common Pattern

Most methods follow the same invocation pattern:

```bash
export PROJECT_ROOT=/path/to/PoseBench

# Hydra-based methods (chai, dynamicbind, rfaa, diffdock, fabind, neuralplexer):
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
data_dir: /home/aoxu/projects/PoseBench/data/runsNposes/
icm_executable: /home/aoxu/icm-3.9-4/icm64
icm_dockscan_path: /home/aoxu/icm-3.9-4/_dockScan
icm_docking_dir: /path/to/PoseBench/forks/ICM
icm_map_dir: /path/to/PoseBench/forks/ICM/ICM_manual_docking_maps/runsNposes
icb_out_dir: /path/to/PoseBench/forks/ICM/inference/runsNposes
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
python3 cogligandbench/models/chai_inference.py dataset=plinder repeat_index=0
```

Key config fields:
```yaml
dataset: plinder
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
**Config:** `configs/model/dynamicbind_inference.yaml` (uses upstream posebench config)

```bash
python3 cogligandbench/models/dynamicbind_inference.py dataset=posebusters_benchmark repeat_index=1
```

Reads protein files from `input_data_dir` and ligand CSV files from `input_ligand_csv_dir`. Pairs proteins and ligands by matching filename stems. Runs DynamicBind via subprocess from `forks/DynamicBind/`.

---

## DiffDock

**Script:** `cogligandbench/models/diffdock_inference.py`
**Config:** `configs/model/diffdock_inference.yaml` (uses upstream posebench config)

```bash
python3 cogligandbench/models/diffdock_inference.py dataset=posebusters_benchmark repeat_index=1
```

Reads from `input_csv_path` (columns: protein path, ligand SMILES). Supports `v1_baseline` and `pocket_only_baseline` flags. Runs DiffDock via subprocess from `forks/DiffDock/`.

---

## FABind

**Script:** `cogligandbench/models/fabind_inference.py`
**Config:** `configs/model/fabind_inference.yaml`

```bash
python3 cogligandbench/models/fabind_inference.py dataset=posebusters_benchmark repeat_index=1
```

---

## NeuralPLexer

**Script:** `cogligandbench/models/neuralplexer_inference.py`
**Config:** `configs/model/neuralplexer_inference.yaml`

```bash
python3 cogligandbench/models/neuralplexer_inference.py dataset=posebusters_benchmark repeat_index=1
```

Supports `no_ilcl` flag to use the rigid docking (`pdbbind_finetuned`) checkpoint instead of the default.

---

## RoseTTAFold-All-Atom (RFAA)

**Script:** `cogligandbench/models/rfaa_inference.py`
**Config:** `configs/model/rfaa_inference.yaml`

```bash
python3 cogligandbench/models/rfaa_inference.py dataset=posebusters_benchmark repeat_index=1
```

---

## Notes

- **Input CSV format** (Vina/GNINA): columns `complex_name`, `protein_path`, `ligand_path`
- **Input directory format** (Chai-1): `<input_dir>/<complex_name>/<complex_name>.fasta`
- **Input directory format** (ICM): `<data_dir>/<complex_name>/<complex_name>.pdb` + `*.sdf`
- **Output poses** are written as ranked SDF files: `<complex_name>_pose<rank>_score<score>.sdf`
- **Timing logs** record `complex_name,elapsed_seconds` per entry
- **FlowDock and AlphaFold 3** do not yet have cogligandbench inference scripts; use their respective `forks/` directories directly

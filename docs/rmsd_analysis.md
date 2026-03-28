# RMSD Analysis Tool for Molecular Docking Results

This tool computes Root Mean Square Deviation (RMSD) between reference ligands and predicted poses from molecular docking methods (Vina and UniDock2), then calculates success rates at different RMSD thresholds.

## Features

- **Multi-method Support**: Analyzes results from both AutoDock Vina and UniDock2
- **Flexible Input Formats**: Handles different output naming conventions automatically
- **Success Rate Calculation**: Computes top-1 and top-5 success rates at customizable RMSD thresholds
- **Comprehensive Reporting**: Generates detailed CSV reports and summary statistics
- **Robust Error Handling**: Continues analysis even if individual systems fail

## Output Formats Supported

### Vina Results
- File pattern: `{system_id}_pose{rank}_score{score:.2f}.sdf`
- Example: `5s9l__1__1.A__1.H_1.I_pose1_score-8.45.sdf`
- Directory structure: `vina_results/{system_id}/pose_files`

### UniDock2 Results  
- File pattern: `rank{rank}.sdf`
- Example: `rank1.sdf`, `rank2.sdf`, etc.
- Directory structure: `unidock2_results/{system_id}/rank_files`

## Installation & Setup

```bash
# Make sure you're in the PoseBench project root
cd /home/aoxu/projects/PoseBench

# The tool uses existing dependencies from the project
# Main requirements: RDKit, pandas, numpy, pyyaml
```

## Usage

### Method 1: Using Configuration File (Recommended)

1. **Edit the configuration file:**
```bash
nano configs/rmsd_analysis.yaml
```

2. **Update paths in the config:**
```yaml
reference_dir: "/path/to/your/reference/dataset"
vina_results: "/path/to/your/vina/results"
unidock2_results: "/path/to/your/unidock2/results"
output_dir: "/path/to/output/directory"
max_poses: 5
thresholds: [2.0, 5.0]
```

3. **Run the analysis:**
```bash
python scripts/rmsd_analysis.py --config configs/rmsd_analysis.yaml
```

### Method 2: Using Command Line Arguments

```bash
python scripts/rmsd_analysis.py \
    --reference_dir /path/to/reference/dataset \
    --vina_results /path/to/vina/results \
    --unidock2_results /path/to/unidock2/results \
    --output_dir /path/to/output \
    --max_poses 5 \
    --thresholds 2.0 5.0
```

### Method 3: Analyze Single Method

```bash
# Analyze only Vina results
python scripts/rmsd_analysis.py \
    --reference_dir /path/to/reference/dataset \
    --vina_results /path/to/vina/results \
    --output_dir /path/to/output

# Analyze only UniDock2 results  
python scripts/rmsd_analysis.py \
    --reference_dir /path/to/reference/dataset \
    --unidock2_results /path/to/unidock2/results \
    --output_dir /path/to/output
```

## Output Files

The tool generates several output files in the specified output directory:

### 1. `detailed_rmsd_results.csv`
Contains RMSD values for every pose analyzed:
```csv
system_id,method,pose_rank,docking_score,rmsd,success_2A,success_5A
5s9l__1__1.A__1.H_1.I,vina,1,-8.45,1.234,True,True
5s9l__1__1.A__1.H_1.I,vina,2,-7.89,2.456,False,True
5s9l__1__1.A__1.H_1.I,unidock2,1,-8.67,0.987,True,True
...
```

### 2. `system_summary.csv`
Per-system summary with best RMSDs and success flags:
```csv
system_id,has_reference,vina_poses,unidock2_poses,vina_best_rmsd,unidock2_best_rmsd,vina_top1_success_2A,vina_top5_success_2A,...
5s9l__1__1.A__1.H_1.I,True,5,5,1.234,0.987,True,True,...
```

### 3. `success_rate_summary.yaml`
Overall success rate statistics:
```yaml
total_systems: 100
vina:
  systems_with_results: 95
  top1_rate_2A: 45.3
  top5_rate_2A: 67.8
  top1_rate_5A: 78.9
  top5_rate_5A: 89.2
unidock2:
  systems_with_results: 98
  top1_rate_2A: 52.1
  top5_rate_5A: 71.4
  ...
```

### 4. `rmsd_analysis.log`
Detailed logging information for debugging and monitoring.

## Understanding the Results

### Success Rate Metrics

- **Top-1 Success Rate**: Percentage of systems where the best-ranking pose (by RMSD) meets the threshold
- **Top-5 Success Rate**: Percentage of systems where any of the top-5 poses (by RMSD) meets the threshold

### Common RMSD Thresholds

- **2.0 Å**: Stringent threshold, indicates very accurate binding pose prediction
- **5.0 Å**: More lenient threshold, often used for initial screening applications

### Interpretation Guidelines

- **RMSD < 2.0 Å**: Excellent prediction, suitable for structure-based drug design
- **RMSD 2.0-5.0 Å**: Good prediction, useful for virtual screening
- **RMSD > 5.0 Å**: Poor prediction, may require method improvement

## Expected Directory Structure

### Reference Dataset
```
reference_dataset/
├── system1_ligand.sdf
├── system2_ligand.sdf
└── ...
```
OR
```
reference_dataset/
├── system1/
│   └── system1_ligand.sdf
├── system2/
│   └── system2_ligand.sdf
└── ...
```

### Vina Results
```
vina_results/
├── system1/
│   ├── system1_pose1_score-8.45.sdf
│   ├── system1_pose2_score-7.89.sdf
│   └── ...
├── system2/
│   └── ...
└── ...
```

### UniDock2 Results
```
unidock2_results/
├── system1/
│   ├── rank1.sdf
│   ├── rank2.sdf
│   └── ...
├── system2/
│   └── ...
└── ...
```

## Testing

Test the tool with sample data:

```bash
python scripts/test_rmsd_analysis.py
```

This will:
1. Check if your data directories are accessible
2. List available systems
3. Test analysis on a single system
4. Provide guidance for running the full analysis

## Troubleshooting

### Common Issues

1. **"No systems found for analysis"**
   - Check that your directory paths are correct
   - Verify the expected file naming conventions
   - Ensure SDF files are readable by RDKit

2. **"Could not load molecule from file"**
   - Check SDF file format and integrity
   - Verify files are not corrupted or empty
   - Try opening files manually with RDKit to test

3. **Import errors**
   - Ensure you're running from the PoseBench project root
   - Check that all dependencies are installed
   - Verify Python path includes the project directory

### Performance Tips

- For large datasets, consider analyzing subsets first
- Monitor memory usage with many systems
- Use `--max_poses` to limit analysis scope
- Check log files for detailed error information

## Advanced Usage

### Custom Thresholds
```bash
python scripts/rmsd_analysis.py \
    --thresholds 1.0 1.5 2.0 3.0 5.0 \
    --config configs/rmsd_analysis.yaml
```

### Analyzing Specific Systems
Modify the reference directory to include only systems of interest, or filter the results CSV files post-analysis.

### Integration with Other Tools
The output CSV files can be easily imported into:
- Pandas for further analysis
- R for statistical analysis
- Excel for manual review
- Visualization tools for plotting

## Example Analysis Workflow

```bash
# 1. Setup directories
mkdir -p /path/to/analysis/output

# 2. Configure analysis
cp configs/rmsd_analysis.yaml my_analysis.yaml
# Edit my_analysis.yaml with your paths

# 3. Test with sample
python scripts/test_rmsd_analysis.py

# 4. Run full analysis
python scripts/rmsd_analysis.py --config my_analysis.yaml

# 5. Review results
ls /path/to/analysis/output/
head /path/to/analysis/output/system_summary.csv
```

This tool provides comprehensive RMSD analysis capabilities for benchmarking molecular docking methods. The flexible input handling and detailed output make it suitable for both research and production use cases.
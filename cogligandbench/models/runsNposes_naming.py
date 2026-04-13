import pandas as pd
import os

inputs_csv = "/home/aoxu/projects/PoseBench/forks/Vina/inference/vina_runsNposes_benchmark_inputs.csv"
data_dir = "/home/aoxu/projects/PoseBench/data/runsNposes"

inputs_df = pd.read_csv(inputs_csv)
for idx, row in inputs_df.iterrows():
    system = row['complex_name']
    selected_ligand_file = row['selected_ligand_file']
    
    system_dir = os.path.join(data_dir, system)
    if os.path.exists(system_dir):
        for file in os.listdir(os.path.join(system_dir, "ligand_files")):
            if file.endswith('.sdf'):
                old_path = os.path.join(system_dir, "ligand_files", file)
                new_path = os.path.join(system_dir, "ligand_files", selected_ligand_file)
                os.rename(old_path, new_path)
                break
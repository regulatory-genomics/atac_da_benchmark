

import yaml
import hashlib
import os

# Function to compute the SHA256 hash of a file
def calculate_sha256(file_path):
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256.update(byte_block)
    return sha256.hexdigest()

# Example list of local data files with metadata
data_files = [
    {
        "file_path": "pbmc_replicate_simulating.h5ad",
        "data_name": "pbmc_replicate_simulating",
        "assay": "simulating"
    }
]

data_files = [
    {
        "file_path": "CD14_Mono_Memory_batch.h5ad",
        "data_name": "CD14_Mono_Memory_batch",
        "assay": "simulating"
    },
    {
        "file_path": "nb_CD14_Mono_Memory.h5ad",
        "data_name": "nb_CD14_Mono_Memory",
        "assay": "simulating"
        
        }
]

data_files = [
    {
        "file_path": "nobatch_CD14_Mono_Memory.h5ad",
        "data_name": "nobatch_CD14_Mono_Memory",
        "assay": "simulating"
    }
]

# Generate data.yaml for each file in the list
for data in data_files:
    file_path = data["file_path"]
    if os.path.exists(file_path):
        data_hash = f"sha256:{calculate_sha256(file_path)}"
        data.update({"data_hash": data_hash})
        data_url = f"{os.path.abspath(file_path)}"
        data.update({"data_url": data_url})

        # Write data.yaml
        yaml_file = f"{os.path.splitext(file_path)[0]}_data.yaml"
        with open(yaml_file, 'w') as file:
            yaml.dump(data_files, file, default_flow_style=False)


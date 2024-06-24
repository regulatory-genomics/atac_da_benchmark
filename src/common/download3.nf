nextflow.enable.dsl=2

include { json_string; json } from './utils.gvy'

process download_dataset {
    container 'jiaqifan/download:v12'
    input:
      val(metadata)
    output:
      tuple val("${json_string(metadata)}"), path("${metadata.data_name}_data.h5ad")

    """
    #!/opt/conda/bin/python
    import shutil
    import os
    import json
    import scanpy as sc
    import pandas as pd
    metadata = json.loads('${json_string(metadata)}')
    file_url = metadata['data_url']
    data=sc.read(file_url)
    data.write(metadata['data_name']+'_data.h5ad')
    """
}

process download_groundtruth {
    container 'jiaqifan/scvi-tools:v8'
    input:
      val(metadata)
    output:
      tuple val("${json_string(metadata)}"), path("data.h5ad") ,path("ground_truth.csv")

    """
    #!/opt/conda/bin/python
    import shutil
    import os
    import json
    import scanpy as sc
    import pandas as pd
    metadata = json.loads('${json_string(metadata)}')
    file_url = metadata['data_url']
    data=sc.read(file_url)
    data.write('data.h5ad')
    data_name = metadata['data_name']
    ground_truth = "/workspace/"+data_name+"_groundTruth.csv"
    truth = pd.read_csv(ground_truth)
    truth.to_csv("ground_truth.csv",index=False)
    """
}

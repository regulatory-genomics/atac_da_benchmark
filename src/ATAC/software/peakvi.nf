nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process da_peakvi {
    container 'jiaqifan/scvi-tools:v3'
    tag "${json(metadata).data_name}"
    cpus 60
    input:
      tuple val(metadata), path("${json(metadata).data_name}_data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'peakvi')}"), path("${json(metadata).data_name}_peakvi_da_region.csv")
    """
    #!/opt/conda/bin/python
    import scvi
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import json
    scvi.settings.seed = 0
    from scipy.sparse import csr_matrix
    data=sc.read(json.loads('${metadata}')['data_name']+"_data.h5ad")
    #min_cells = int(data.shape[0] * 0.05)
    #sc.pp.filter_genes(data , min_cells=min_cells)
    if 'batch' in data.obs:
        #   batch_key = data.obs['batch']
        scvi.model.PEAKVI.setup_anndata(data, batch_key='batch')
    else:
        scvi.model.PEAKVI.setup_anndata(data)
    
    pvi = scvi.model.PEAKVI(data, n_latent=30)
    pvi.train()
    celltype=data.obs['cell_type'].unique()
    group=data.obs['compare'].unique()
    g1='ref'
    g2='case'
    index = 0
    da_list={}
    while index < len(celltype):
         g1_index = (data.obs["cell_type"] == celltype[index]) & (data.obs["compare"] == g1)
         g2_index = (data.obs["cell_type"] == celltype[index]) & (data.obs["compare"] == g2)
         da_res11 = pvi.differential_accessibility(
    	        	idx1=g1_index,  idx2=g2_index, 
                two_sided=True,batch_correction=True
		      )
         da_list[celltype[index]]=da_res11
         index = index+1

    combined_dars = []

    for cell_type, df in da_list.items():
      df["cell_type"] = cell_type
      combined_dars.append(df) 
    # Store the matrix as a layer
    combined_dars_df = pd.concat(combined_dars)
    combined_dars_df.reset_index(inplace=True)

    #adata = ad.AnnData(
    #    X=np.empty((len(combined_dars_df), 0)),  # Empty matrix
    #    obs=combined_dars_df
    #)
    # Save to .h5ad format
    #adata.write("gene_cell_matrices.h5ad")
    combined_dars_df.to_csv(json.loads('${metadata}')['data_name']+"_peakvi_da_region.csv", index=False)



    """
}

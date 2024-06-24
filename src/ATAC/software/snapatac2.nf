nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process da_snapatac2 {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    cpus 60
    input:
      tuple val(metadata), path("${json(metadata).data_name}_data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'snapatac2')}"), path("${json(metadata).data_name}_snapatac2_da_region.csv")
    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import polars as pl
    import json
    data=sc.read(json.loads('${metadata}')['data_name']+"_data.h5ad")
    #min_cells = int(data.shape[0] * 0.05)
    #sc.pp.filter_genes(data , min_cells=min_cells)
    
    celltype=data.obs['cell_type'].unique()
    g1='ref'
    g2='case'
    index = 0
    da_list={}
    while index < len(celltype):
         g1_index = (data.obs["cell_type"] == celltype[index]) & (data.obs["compare"] == g1)
         g2_index = (data.obs["cell_type"] == celltype[index]) & (data.obs["compare"] == g2)
         if len(data.obs['batch'].unique()) > 1:
          da_res11 = snap.tl.diff_test( data,
                cell_group1 = g1_index,  cell_group2 = g2_index, 
                  min_log_fc=0.0005, min_pct=0.0001
            )
         else:
          da_res11 = snap.tl.diff_test( data,
                cell_group1 = g1_index,  cell_group2 = g2_index, 
                  min_log_fc=0.0005, min_pct=0.0001
            )
         da_list[celltype[index]]=da_res11
         index = index+1

    combined_dars = []

    for cell_type, df in da_list.items():
    	modified_df = df.clone()
    	modified_df = modified_df.with_columns(pl.lit(cell_type).alias("cell_type"))
    	combined_dars.append(modified_df) 
    combined_dars_t = pl.concat(combined_dars)
    combined_dars_t = combined_dars_t.rename({'feature name': 'index', 'p-value': 'pvalue', 'adjusted p-value': 'padj'})
    combined_dars_t.write_csv(json.loads('${metadata}')['data_name']+"_snapatac2_da_region.csv")



    """
}

nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'
process da_archr {
    container 'jiaqifan/archr:v5'
    tag "${json(metadata).data_name}"
    cpus 60

    when: is_included("scada", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("${json(metadata).data_name}_data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'archr')}"), path("${json(metadata).data_name}_archr_da_region.csv")

    """
    #!/usr/local/bin/Rscript
    library("rhdf5")
    library(jsonlite) 
    library(Matrix)
    library("presto")
    library('S4Vectors')
    library("nabor")
    json_data <- jsonlite::fromJSON('${metadata}')
    file_name <- json_data\$data_name
    read_h5ad_data <- function(h5ad_file, group = "X") {
      file <- H5Fopen(h5ad_file, flags = "H5F_ACC_RDONLY")
      x <- H5Gopen(file, group)

      indices <- H5Dread(H5Dopen(x, "indices"))
      indptr <- H5Dread(H5Dopen(x, "indptr"))
      data <- H5Dread(H5Dopen(x, "data"))
      shape <- rev(H5Aread(H5Aopen(x, "shape")))
      sparse_matrix <- sparseMatrix(i = indices+1, p = indptr, x = as.numeric(data), index1 = TRUE, repr = "C")
      H5Gclose(x)
      obs_names <- H5Dread(H5Dopen(file, "obs/_index"))
      var_names <- H5Dread(H5Dopen(file, "var/_index"))
      H5Fclose(file)
      if(length(var_names)==ncol(sparse_matrix)){
      rownames(sparse_matrix) <- obs_names
      colnames(sparse_matrix) <- var_names}else{
      colnames(sparse_matrix) <- obs_names
      rownames(sparse_matrix) <- var_names
    }
  
      return(sparse_matrix)
    }

    h5ad_file=paste(file_name,'data.h5ad',sep="_")
    data_matrix <- read_h5ad_data(h5ad_file)


    read_h5ad_obs <- function(h5ad_file, group_t = "obs", columns = NULL) {
      file <- H5Fopen(h5ad_file, flags = "H5F_ACC_RDONLY")
  
      datasets <- h5ls(file, recursive = TRUE)
      obs_datasets <- datasets[datasets\$group == paste0("/", group_t), ]

      if (is.null(columns)) {
        columns <- obs_datasets
      } else {
        columns <- intersect(columns, obs_datasets\$name)
      }
  
      obs_index <- H5Dread(H5Dopen(file, paste0(group_t, "/_index")))
  
      obs_data <- lapply(columns, function(col) {
      if (obs_datasets[which(obs_datasets\$name==col),"dim"]==""){
        group_codes <- H5Dread(H5Dopen(file, paste0(group_t,"/",col,"/codes")))
        group_categories <- H5Dread(H5Dopen(file, paste0(group_t,"/",col,"/categories")))
        group_mapping <- setNames(group_categories, seq_along(group_categories) - 1)
        group <- group_mapping[as.character(group_codes)]}else{
      group <- H5Dread(H5Dopen(file, paste0(group_t,"/",col)))

      }
  
      })
  
      obs_df <- as.data.frame(obs_data, stringsAsFactors = FALSE)
      colnames(obs_df) <- columns
      rownames(obs_df) <- obs_index
      H5Fclose(file)
  
      return(obs_df)
    }

    obs_df <- read_h5ad_obs(h5ad_file,columns=c("batch","cell_type","compare"))
    source('/other_methods_for_differential_updated.R') 
    g1='ref'
    g2='case'
    res_sum<-NULL
    if(nrow(obs_df)==ncol(data_matrix)){
            data_matrix=t(data_matrix)}
    cell_type=unique(obs_df\$cell_type)
    for ( type in cell_type){
	      mat_g1=data_matrix[which(obs_df\$compare==g1 & obs_df\$cell_type==type),]
        mat_g2=data_matrix[which(obs_df\$compare==g2 & obs_df\$cell_type==type),]
        res = archR_method(t(mat_g1),t(mat_g2))
	      res_d = data.frame(pval=res,index=colnames(mat_g1),cell_type=type)
        res_sum<-rbind(res_sum,res_d)}
    write.table(res_sum,paste(file_name,"_archr_da_region.csv",sep=""),col.names=T,sep=",",quote=F,row.names=F)
     """
}


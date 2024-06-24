nextflow.enable.dsl = 2

include { is_included; add_meta; json; } from '../../common/utils.gvy'
process da_deseq2 {
    container 'jiaqifan/deseq2:v1'
    tag "${json(metadata).data_name}"
    cpus 60

    when: is_included('scada', params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("${json(metadata).data_name}_data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'deseq2')}"), path("${json(metadata).data_name}_deseq2_da_region.csv")

    """
    #!/usr/local/bin/Rscript
    library("rhdf5")
    library(Matrix)
    library("DESeq2")
    library(jsonlite)
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
    ##which(colSums(data_matrix)==0)
    if(nrow(obs_df)==ncol(data_matrix)){
            data_matrix=t(data_matrix)}
    cell_type_list=unique(obs_df\$cell_type)
    merge_sample=function(matrix,sample_index){
      bulk=NULL
      for(i in 1:(length(sample_index)-1)){
      pseudo=colSums(matrix[sample_index[i]:sample_index[i+1],])
      bulk=cbind(bulk,pseudo)
      }
      bulk
    }
    res_sum<-NULL
    for(type in cell_type_list){
      obs_df_ref=obs_df[intersect(which(obs_df\$cell_type==type),which(obs_df\$compare=="ref")),]
      obs_df_case=obs_df[intersect(which(obs_df\$cell_type==type),which(obs_df\$compare=="case")),]
      ref_n=round(nrow(obs_df_ref)/4)
      ref_sample=seq(1,nrow(obs_df_ref),ref_n)
      ref_bulk=merge_sample(data_matrix[rownames(obs_df_ref),],ref_sample)
      colnames(ref_bulk)=paste("ref",1:ncol(ref_bulk),sep="_")

      case_n=round(nrow(obs_df_case)/4)
      case_sample=seq(1,nrow(obs_df_case),case_n)
      case_bulk=merge_sample(data_matrix[rownames(obs_df_case),],case_sample)
      colnames(case_bulk)=paste("case",1:ncol(case_bulk),sep="_")
      coldata=data.frame(cell_type=type,compare=c(rep("ref",ncol(ref_bulk)),
      rep("case",ncol(case_bulk))),batch=c(1:ncol(ref_bulk),1:ncol(case_bulk))
      )
      rownames(coldata)=c(colnames(ref_bulk),colnames(case_bulk))

      dds <- DESeqDataSetFromMatrix(countData = cbind(ref_bulk,case_bulk) ,
                                colData = coldata,
                                design= ~ compare)

      dds <- DESeq(dds)
      #resultsNames(dds) # lists the coefficients
      res <- results(dds, contrast=c("compare", "case", "ref"))
      dm <- data.frame(res)
      dm\$cell_type=type
      dm\$index=rownames(dm)
      res_sum<-rbind(res_sum,dm)
    }
    write.table(res_sum,paste(file_name,"_deseq2_da_region.csv",sep=""),col.names=T,sep=",",quote=F,row.names=F)
     """
}

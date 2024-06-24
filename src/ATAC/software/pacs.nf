nextflow.enable.dsl=2
include { is_included; add_meta; json; } from '../../common/utils.gvy'
process da_pacs {
    container 'jiaqifan/pacs:v3'
    tag "${json(metadata).data_name}"
    cpus 60

    when: is_included("pacs", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("${json(metadata).data_name}_data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'pacs')}"), path("${json(metadata).data_name}_pacs_da_region.csv")

    """

    #!/usr/bin/env Rscript
    library(PACS)
    library(Matrix)
    library(PICsnATAC)
    library(rhdf5)
    library(jsonlite)
    json_data <- jsonlite::fromJSON('${metadata}')
    file_name <- json_data\$data_name
    #BiocManager::install("rhdf5")

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

    h5ad_file=paste(file_name,'data.h5ad',sep="_")
    data_matrix <- read_h5ad_data(h5ad_file)
    obs_df <- read_h5ad_obs(h5ad_file,columns=c("batch","cell_type","compare"))
    #mat=readMM('GSE129785_scATAC-Hematopoiesis-All.mtx.gz')
    #metadata=read.table("GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt",header=T)
    if(nrow(obs_df)==ncol(data_matrix)){
            data_matrix=t(data_matrix)}

    r_by_ct_out <- get_r_by_ct_mat_pq(
	    cell_type_set =unique(obs_df\$cell_type) ,
  	  r_by_c = t(data_matrix),
 	    cell_type_labels = obs_df\$cell_type,
  	  n_features_per_cell = dim(t(data_matrix))[1]
    )
    data_matrix_d=as.matrix(data_matrix)
    print(1)
    #######################-----------##
    ## wrapper for covariate matrix

    pacs_test_sparse <- function(covariate_meta.data, formula_full,
                             formula_null, pic_matrix,
                             n_peaks_per_round = NULL,
                             T_proportion_cutoff = 0.2,
                             cap_rates, par_initial_null = NULL,
                             par_initial_full = NULL, n_cores = 1,
                             verbose = TRUE) {
  	## check the data
 	  n_cell <- ncol(pic_matrix)
	  n_peaks <- nrow(pic_matrix)
  	p_names <- rownames(pic_matrix)

  	if (nrow(covariate_meta.data) != n_cell) {
    		stop("number of cells do not match between meta.data and data matrix")
 	  }

	  if (length(cap_rates) != n_cell) {
    		stop("number of cells do not match between cap_rates and data matrix")
 	  }

    if (is.null(rownames(pic_matrix))) {
          message("peak names not supplied, set to f_1 to f_n")
          rownames(pic_matrix) <- paste("f", 1:n_peaks, sep = "_")
    }

  	if (is.null(n_peaks_per_round)) {
    		n_peaks_per_round <- min(floor(2^29 / n_cell), n_peaks) ## to be safer
  	}


  	if (inherits(pic_matrix, "Matrix")) {
    		## counts that are >= 2
    		pic_matrix_2 <- pic_matrix
    		pic_matrix_2@x[pic_matrix_2@x == 1] <- 0

    		pic_matrix_2 <- drop0(pic_matrix_2)
    		pic_matrix_2@x <- rep(1, length = length(pic_matrix_2@x))

    		## binarized matrix
    		pic_matrixbin <- pic_matrix
    		pic_matrixbin@x <- rep(1, length = length(pic_matrixbin@x))
  	} else {
    		## counts that are >= 2
    		pic_matrix_2 <- Matrix::Matrix(pic_matrix, sparse = TRUE)
    		pic_matrix_2@x[pic_matrix_2@x == 1] <- 0

    		pic_matrix_2 <- Matrix::drop0(pic_matrix_2)
    		pic_matrix_2@x <- rep(1, length = length(pic_matrix_2@x))

    		## binarized matrix
    		pic_matrixbin <- Matrix::Matrix(pic_matrix, sparse = TRUE)
    		pic_matrixbin@x <- rep(1, length = length(pic_matrixbin@x))
    		pic_matrixbin = as.matrix(pic_matrixbin)
  	}

    ## calculate the proportion of counts >= 2
    rs <- rowSums(pic_matrixbin)
    rs2 <- rowSums(pic_matrix_2)

    p_2 <- rs2 / rs ## proportion of counts >= 2
    p_2[is.na(p_2)] <- 0

    n_p_2 <- sum(p_2 >= T_proportion_cutoff)
    n_p_b <- sum(p_2 < T_proportion_cutoff)

    if (verbose) {
      print(paste(n_p_2, "peaks consider cumulative logit models", sep = " "))
      print(paste(n_p_b, "peaks consider logit models", sep = " "))
    }

    rm(pic_matrix_2)
    gc(verbose = FALSE)

    ## select features based on the proportion of 2
    f_sel <- names(p_2)[p_2 >= T_proportion_cutoff]
    f_b_sel <- names(p_2)[p_2 < T_proportion_cutoff]

    if (n_p_2 >= 1) {
      pic_matrix <- pic_matrix[f_sel, , drop = FALSE]

      n_iters <- ceiling(n_p_2 / n_peaks_per_round)
      p_cumu <- rep(list(), length = n_iters)

      for (jj in 1:n_iters) {
        peak_start <- (jj - 1) * n_peaks_per_round + 1
        peak_end <- min(n_p_2, jj * n_peaks_per_round)

        pic_dense <- as.matrix(pic_matrix[peak_start:peak_end, ])

        print("pacs computing cumulative logit part")

        p_cumu[[jj]] <- pacs_test_cumu(
          covariate_meta.data = covariate_meta.data, max_T = 2,
          formula_full = formula_full,
          formula_null = formula_null,
          pic_matrix = pic_dense,
          cap_rates = cap_rates, n_cores = n_cores,
          par_initial_null = par_initial_null,
          par_initial_full = par_initial_full
        )
      # toc()
      }
      rm(pic_dense)
    }
    ## we do not need the quantitative part of the matrix any more
    rm(pic_matrix)
    gc(verbose = FALSE)


    ## compute the binary part
    if (n_p_b >= 1) {
      pic_matrixbin <- pic_matrixbin[f_b_sel, , drop = FALSE]

      n_iters_b <- ceiling(n_p_b / n_peaks_per_round)
      p_logit <- rep(list(), length = n_iters_b)

      for (jj in 1:n_iters_b) {
        peak_start <- (jj - 1) * n_peaks_per_round + 1
        peak_end <- min(n_p_b, jj * n_peaks_per_round)

        pic_dense <- as.matrix(pic_matrixbin[peak_start:peak_end, ,drop = FALSE ])

        print("pacs computing logit part")

        p_logit[[jj]] <- pacs_test_logit(
          covariate_meta.data = covariate_meta.data,
          formula_full = formula_full,
          formula_null = formula_null,
          pic_matrix = pic_dense,
          cap_rates = cap_rates, n_cores = n_cores,
          par_initial_null = par_initial_null,
          par_initial_full = par_initial_full
        )
        #toc()
      }
      rm(pic_dense)
      gc(verbose = FALSE)
    }

    ## organize the p values as well as convergence status
    if (n_p_2 >= 1 && n_p_b >= 1) {
      ## cumulative part
      if (n_iters > 1) {
        ## p values
        p_cumu_list <- sapply(p_cumu, function(x) x\$pacs_p_val)
        p_val_cumu <- unlist(p_cumu_list)

        ## convergence status
        conv_cumu_list <- sapply(p_cumu, function(x) {
          matrix(x\$pacs_converged, ncol = 2)
        })
        conv_cumu <- do.call(rbind, conv_cumu_list)
      } else {
        ## p values
        p_cumu_list <- sapply(p_cumu, function(x) x\$pacs_p_val)
        p_val_cumu <- p_cumu_list[, 1, drop = TRUE]

        ## convergence status
        conv_cumu_list <- sapply(p_cumu, function(x) x\$pacs_converged)
        conv_cumu <- matrix(conv_cumu_list, ncol = 2)
        rownames(conv_cumu) <- names(p_val_cumu)
      }

      ## logit part
      if (n_iters_b > 1) {
        ## p values
        p_logit_list <- sapply(p_logit, function(x) x\$pacs_p_val)
        p_val_logit <- unlist(p_logit_list)

        ## convergence status
        conv_logit_list <-
          sapply(p_logit, function(x) matrix(x\$pacs_converged, ncol = 2))
        conver_logit <- do.call(rbind, conv_logit_list)
      } else {
        ## p values
        p_logit_list <- sapply(p_logit, function(x) x\$pacs_p_val)
        p_val_logit <- p_logit_list[, 1, drop = TRUE]

        ## convergence status
        conv_logit_list <- sapply(p_logit, function(x) x\$pacs_converged)
        conver_logit <- matrix(conv_logit_list, ncol = 2)
        rownames(conver_logit) <- names(p_val_logit)
      }

      p_val <- c(p_val_cumu, p_val_logit)[p_names]
      convergence <- rbind(conv_cumu, conver_logit)[p_names, ]
    } else if (n_p_2 == 0) {
        ## p values
        p_logit_list <- sapply(p_logit, function(x) x\$pacs_p_val)
        p_val_logit <- unlist(p_logit_list)
        p_val <- p_val_logit[p_names]

        ## convergence status
        conv_logit_list <-
          sapply(p_logit, function(x) matrix(x\$pacs_converged, ncol = 2))
        conver_logit <- do.call(rbind, conv_logit_list)
        rownames(conver_logit) <- names(p_val_logit)
        convergence <- conver_logit[p_names, ]

    } else if (n_p_b == 0) {
        ## p values
        p_cumu_list <- sapply(p_cumu, function(x) x\$pacs_p_val)
        p_val_cumu <- unlist(p_cumu_list)
        p_val <- p_val_cumu[p_names]

        ## convergence status
        conv_cumu_list <-
          sapply(p_cumu, function(x) matrix(x\$pacs_converged, ncol = 2))
        conv_cumu <- do.call(rbind, conv_cumu_list)
        rownames(conver_logit) <- names(p_val_logit)
        convergence <- conv_cumu[p_names, ]
    }


    return(list(pacs_converged = convergence, pacs_p_val = p_val))
  }

  ####
  res_sum<-NULL
  print(2)
  celltype_list=unique(obs_df\$cell_type)
  for (i in celltype_list){
  obs_df_type=obs_df[which(obs_df\$cell_type==i),]
  data_matrix_type=data_matrix_d[rownames(obs_df_type),]
  rate_type=r_by_ct_out\$q_vec[rownames(obs_df_type)]
  if(length(unique(obs_df_type\$batch))>1){
  p_vals <- pacs_test_sparse(
          covariate_meta.data = obs_df_type,
          formula_full = ~ factor(compare)+factor(batch),
          formula_null = ~ factor(batch),
          pic_matrix = t(data_matrix_type),
          cap_rates = rate_type)

  }else{
      p_vals <- pacs_test_sparse(
    covariate_meta.data = obs_df_type,
    formula_full = ~ factor(compare),
    formula_null = ~ 1,
    pic_matrix = t(data_matrix_type),
    cap_rates = rate_type
    )}
  res=data.frame(index=names(p_vals\$pacs_p_val),pval=p_vals\$pacs_p_val,cell_type=i)
  res_sum<-rbind(res_sum,res)

  }
  write.table(res_sum,paste(file_name,"_pacs_da_region.csv",sep=""),col.names=T,sep=",",row.names=F,quote=F)

    """
}

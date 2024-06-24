nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'
process da_scada {
    container 'jiaqifan/scada:v4'
    tag "${json(metadata).data_name}"
    cpus 60

    when: is_included("scada", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("${json(metadata).data_name}_data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'scada')}"), path("${json(metadata).data_name}_scada_da_region.csv")

    """
    #!/usr/local/bin/Rscript
    library(scaDA)
    library(parallel)
   
    library("rhdf5")
    library(jsonlite) 
    library(Matrix)
    json_data <- jsonlite::fromJSON('${metadata}')
    file_name <- json_data\$data_name
    zinb.loglink <- function(counts,p,u,k){
      counts <- as.numeric(counts)
      dens <- numeric(length(counts))
      for (i in 1:length(counts)){
        if (counts[i]==0){
          dens[i] <- p+(1-p)*(1/(1+k*u))^(1/k)} else {
          g <- gamma(counts[i]+1/k)/(gamma(1/k)*gamma(counts[i]+1))
          dens[i] <- (1-p)*g*(1/(1+k*u))^(1/k)*((k*u)/(1+k*u))^counts[i]
          }
        dens <- unlist(dens)
      } 
      loglink <- sum(log(dens))
      return(loglink)
      }


    optParamsParallel <- function(object) {
    message("start optimize parameter estimates")

    count <- object@count
    group.1.loc <- object@params\$g1
    group.2.loc <- object@params\$g2
    dat <- count[,c(group.1.loc, group.2.loc)]
    npeak <- dim(dat)[1]
    nsam <- dim(dat)[2]
    nsam1 <- length(group.1.loc)
    nsam2 <- length(group.2.loc)
    poolCol <- c(1:nsam)
    cond1Col <- c(1:nsam1)
    cond2Col <- c((nsam1+1):nsam)

    est_params_cell1 <- data.frame(object@params\$param_g1)
    est_params_cell2 <- data.frame(object@params\$param_g2)

    tol <- 1e-3
    nitr <- 10

    #  no_cores <- parallel::detectCores() - 1
    no_cores <- 60
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterExport(cl, c("dat", "cond1Col", "cond2Col", "nitr", "tol", "est_params_cell1", "est_params_cell2", "zinb.loglink"),envir = environment())

    results_c1 <- parallel::parLapply(cl, 1:npeak, function(i) {
      counts <- dat[i, cond1Col]
      prev <- est_params_cell1[i,]\$p0
      nb_mu <- est_params_cell1[i,]\$mu
      max.mu <- max(est_params_cell1\$mu)
      nb_phi <- est_params_cell1[i,]\$phi
      for (k in 1:nitr) {
        prev0 <- prev
        nb_mu0 <- nb_mu
        nb_mu <- optimise(zinb.loglink, c(0.01, max.mu), tol = 1e-4, maximum = TRUE, counts = counts, p = prev0, k = nb_phi)\$maximum
        prev <- optimise(zinb.loglink, c(0.01, 1), tol = 1e-4, maximum = TRUE, counts = counts, u = nb_mu, k = nb_phi)\$maximum
        if (abs(nb_mu0 - nb_mu) / abs(nb_mu0) < tol && abs(prev0 - prev) / abs(prev0) < tol) {
          break
        }
        }
       return(list(mu = nb_mu, prev = prev))
      })
      # Process results for condition 1
      mu_opt_c1 <- sapply(results_c1, function(x) x\$mu)
      prev_opt_c1 <- sapply(results_c1, function(x) x\$prev)

      # Parallelize optimization for condition 2
        results_c2 <- parallel::parLapply(cl, 1:npeak, function(i) {
        counts <- dat[i, cond2Col]
        prev <- est_params_cell2[i,]\$p0
        nb_mu <- est_params_cell2[i,]\$mu
        max.mu <- max(est_params_cell2\$mu)
        nb_phi <- est_params_cell2[i,]\$phi
        for (k in 1:nitr) {
          prev0 <- prev
          nb_mu0 <- nb_mu
          nb_mu <- optimise(zinb.loglink, c(0.01, max.mu), tol = 1e-4, maximum = TRUE, counts = counts, p = prev0, k = nb_phi)\$maximum
          prev <- optimise(zinb.loglink, c(0.01, 1), tol = 1e-4, maximum = TRUE, counts = counts, u = nb_mu, k = nb_phi)\$maximum
          if (abs(nb_mu0 - nb_mu) / abs(nb_mu0) < tol && abs(prev0 - prev) / abs(prev0) < tol) {
            break
          }
        }
        return(list(mu = nb_mu, prev = prev))
        })

      # Process results for condition 2
        mu_opt_c2 <- sapply(results_c2, function(x) x\$mu)
        prev_opt_c2 <- sapply(results_c2, function(x) x\$prev)

        parallel::stopCluster(cl)

      # updates parameter estimates
        est_params_cell1\$mu <- mu_opt_c1
        est_params_cell1\$p0 <- prev_opt_c1
        est_params_cell2\$mu <- mu_opt_c2
        est_params_cell2\$p0 <- prev_opt_c2
        est_params_pooled <- object@params\$param_pooled

      ## testing using shrinked phi and optimized mu and prev
        pval_zinb_shrink_opt <- NULL
        tstats <- NULL
        for (i in 1:npeak){
          counts <- dat[i,poolCol]
          prev <- est_params_pooled[i,]\$p0
          nb_mu <- est_params_pooled[i,]\$mu
          nb_phi_shrink <- est_params_pooled[i,]\$phi
          logL_null <- zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)

          counts <- dat[i,cond1Col]
          prev <- est_params_cell1[i,]\$p0
          nb_mu <- est_params_cell1[i,]\$mu
          nb_phi_shrink <- est_params_cell1[i,]\$phi
          logL_alter_1 <- zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)

          counts <- dat[i,cond2Col]
          prev <- est_params_cell2[i,]\$p0
          nb_mu <- est_params_cell2[i,]\$mu
          nb_phi_shrink <- est_params_cell2[i,]\$phi
          logL_alter_2 <- zinb.loglink(counts=counts,p=prev,u=nb_mu,k=nb_phi_shrink)
          logL_alter <- logL_alter_1+logL_alter_2

          test.stats <- -2*(logL_null - logL_alter)

          pvl <- pchisq(test.stats, df=3, lower.tail = FALSE)
          pval_zinb_shrink_opt <- c(pval_zinb_shrink_opt,pvl)
          tstats <- c(tstats,test.stats)
          }
          result <- data.frame(tstats=tstats, pval=pval_zinb_shrink_opt)

          result\$FDR <- p.adjust(result\$pval,method='fdr')
      # calculate fold change
          m1 <- (1-est_params_cell1\$p0)*est_params_cell1\$mu
          m2 <- (1-est_params_cell2\$p0)*est_params_cell2\$mu
          foch <- m2/m1
          result\$log2fc <- log(foch,2)
      # update params estimates
          object@params <- list(g1 = group.1.loc,
                            g2 = group.2.loc,
                            param_g1 = est_params_cell1,
                            param_g2 = est_params_cell2)
      # save DA test results to object@result slot
        object@result <- result
        return(object)
        }

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

    h5ad_file <- paste(file_name,"data.h5ad",sep="_")
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
    if(nrow(obs_df)==ncol(data_matrix)){
            data_matrix=t(data_matrix)}
    g1='ref'
    g2='case'
    res_sum<-NULL
    cell_type=unique(obs_df\$cell_type)
    for ( type in cell_type){
	    mat_g1=data_matrix[which(obs_df\$compare==g1 & obs_df\$cell_type==type),]
      mat_g2=data_matrix[which(obs_df\$compare==g2 & obs_df\$cell_type==type),]    	
	    counts_sum_suni=apply(mat_g1,2,function(x){length(which(x>0))})
    	counts_sum_isoc=apply(mat_g2,2,function(x){length(which(x>0))})
    	plot_region=intersect(which(counts_sum_suni>0),which(counts_sum_isoc>0))
	    mat_filterregion=data_matrix[c(rownames(mat_g1),rownames(mat_g2)),plot_region]
   	  #mat_filterregion=data_matrix[c(rownames(mat_g1),rownames(mat_g2)),]
      obs_df_sub=obs_df[rownames(mat_filterregion),]
	    scaDA.obj <- scaDAdatasetFromMatrix(count = t(as.matrix(mat_filterregion)), colData = data.frame(obs_df_sub\$compare))
	    scaDA.obj <- estParams(scaDA.obj, group.1 = g1, group.2 = g2)
    	print(4)
    	scaDA.obj <- shrinkDisp(scaDA.obj)
      # scaDA.obj <- optParams(scaDA.obj)
    	print(5)
    	scaDA.obj <- optParamsParallel(scaDA.obj)
    	print(6)
    	results3 = scaDA.obj@result
	    results3\$cell_type=type
      results3\$index=colnames(mat_filterregion)
      res_sum<-rbind(res_sum,results3)}
      write.table(res_sum,paste(file_name,"_scada_da_region.csv",sep=""),col.names=T,sep=",",quote=F,row.names=F)
    """
}


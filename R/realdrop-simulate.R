
#' @export
#' 

RealDropSimulate <- function(params = newRealDropParams(), realdata,
                          verbose = TRUE, ...){
  
  checkmate::assertClass(params, "RealDropParams")
  
  if (verbose) {message("Getting parameters...")}
  params <- setParams(params, ...)
  params <- expandParams(params)
  validObject(params)
  
  # Set random seed
  seed <- getParam(params, "seed")
  set.seed(seed)
  
  # Get the parameters we are going to use
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  
  genethred <- getParam(params, "genethred")
  cellthred <- getParam(params, "cellthred")

  cellthr = quantile(colSums(realdata), cellthred)
  genethr = quantile(rowSums(realdata>0), genethred)
  cellindex = which(colSums(realdata)>cellthr)
  geneindex = which(rowSums(realdata>0)>genethr)
  
  if(nGenes > length(geneindex)){
    cat("Warning: the requested nGenes are too much with regard gene threshod", genethred,", modifying to ", length(geneindex),"... \n")
    nGenes = length(geneindex)
    setParams(params, nGenes = nGenes)
  }
  if(nCells > length(cellindex)){
    cat("Warning: the requested nCells are too much with regard cell threshod", cellthred,", modifying to ", length(cellindex),"... \n")
    nCells = length(cellindex)
    setParams(params, nCells = nCells)
  }
  
  if (verbose) {message("Creating simulation object...")}
  
  true_data = realdata[sample(geneindex, nGenes), sample(cellindex, nCells)]

  # Set up name vectors
  cell.names <- colnames(true_data)
  gene.names <- rownames(true_data)
  
  # Create SingleCellExperiment to store simulation
  cells <-  data.frame(Cell = cell.names)
  rownames(cells) <- cell.names
  features <- data.frame(Gene = gene.names)
  rownames(features) <- gene.names
  sim <- SingleCellExperiment(rowData = features, colData = cells,
                              metadata = list(Params = params))
  
  assays(sim)$TrueCounts = true_data
  
  trials <- getParam(params, "trials")
  dirname <- getParam(params, "dirname")
  
  if(trials>1){
    if(!dir.exists(dirname)){
      message("Warning: detect calls for multiple trials, but no directory to save files...")
      message("Modify: only do the simulations with one trial...")
      params <- setParams(params, trials = 1)
      metadata(sim)$Params= params
    }
  }
  
  if(!dir.exists(dirname)){
    if (verbose) message("Adding technical noise ...")
    sim<-ObservedCounts(sim, "") 
    sim<-RealSimDropout(sim, "")
    return(sim)
  }
  else{
    if (verbose) message("Adding technical noise ...")
    trials<-getParam(metadata(sim)$Params, "trials")
    if(trials>1){
      message("Starting multiple trials.....")
      for(trial in 1:trials){
        progress(trial, progress.bar = TRUE)
        Sys.sleep(0.01)
        sim<-ObservedCounts(sim, trial) 
        sim<-RealSimDropout(sim, trial)
      }
    }
    else{
      sim<-ObservedCounts(sim, "") 
      sim<-RealSimDropout(sim, "")
    }
    if(dir.exists(dirname))saveRDS(sim, paste0(dirname, "sim.rds"))
    return(sim)
  }
  
}



#' @export

RealSimDropout<-function(sim, trial){
  
  params <- metadata(sim)$Params
  down.meanvec <- getParam(params, "down.mean")
  down.rate <- getParam(params, "down.rate")
  
  dropout.midvec <- getParam(params, "dropout.mid")
  dropout.shape <- getParam(params, "dropout.shape")
  
  dirname <- getParam(params, "dirname")
  
  true_data = assays(sim)$TrueCounts
  nc = ncol(true_data)
  ng = nrow(true_data)
  
  if((length(down.meanvec)>1 | length(dropout.midvec)>1) & !dir.exists(dirname)){
    message("Warning: detect calls for multiple configurations, but no directory to save files...")
    message("Modify: only do the simulations with one configuration...")
    down.meanvec = down.meanvec[1]
    dropout.midvec = dropout.midvec[1]
    setParams(params, dropout.mid = dropout.midvec, down.mean = down.meanvec)
  }
  
  for(down.mean in down.meanvec){
    tau = rgamma(nc, down.mean, down.rate)
    down_data = round(t(t(true_data)*tau), 0)
    
    for(dropout.mid in dropout.midvec){
      drop.prob <- sapply(1:nc, function(idx) {
          eta <- log(down_data[, idx]+1)
          return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
        })
      keep.prob <- 1 - drop.prob
      keep <- matrix(rbinom(nc*ng, 1, keep.prob), nrow = ng, ncol = nc)
      raw_data <- down_data * keep
      raw_data[is.na(raw_data)] = 0 
      rownames(raw_data) = rownames(true_data)
      colnames(raw_data) = colnames(true_data)
      if(dir.exists(dirname))saveRDS(raw_data, paste0(dirname, trial,"_downme", down.mean, "_dm", dm, "_Raw.rds"))
    }
  }
  assays(sim)$counts = raw_data
  metadata(sim)$Params = params
  return(sim)
}


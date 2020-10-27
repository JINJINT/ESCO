#' ESCO simulation
#'
#' Simulate count data from a fictional single-cell RNA-seq experiment using
#' the hierarchical copula model.
#'
#' @param type which type of heterogenounity to use. Options are :
#'        "single" which produces a single population;
#'        "group" which produces distinct groups;
#'        "tree"  which produces distinct groups but admits a tree structure;
#'        "traj"  which produces distinct groups but admits a smooth trajectory structure.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{\link{params}}.
#'
#' @details
#' Parameters can be set in a variety of ways. If no parameters are provided
#' the default parameters are used. We adopts the way in splatter to set up parameters, particularly, 
#' any parameters in \code{\link{params}} can be
#' overridden by supplying additional arguments through a call to
#' \code{\link{setParams}}. This design allows the user flexibility in
#' how they supply parameters and allows small adjustments without creating a
#' new \code{escoParams} object. See examples for a demonstration of how this
#' can be used.
#'
#' The simulation involves the following steps:
#' \enumerate{
#'     \item Set up simulation object
#'     \item Simulate library sizes
#'     \item Simulate base gene means
#'     \item Simulate groups / tree / path structure differential factors
#'     \item Simulate mean variance relationship adjusted cell means
#'     \item Simulate true counts with / without correlation
#'     \item Simulate techinical noise (zeroinflation / downsampling)
#'     \item Create final dataset
#' }
#'
#' The final output is a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object that
#' contains the simulated counts but also the values for various intermediate
#' steps, and the paramters used to generate the simulation.
#' These are stored in the \code{\link{colData}} (for cell specific
#' information), \code{\link{rowData}} (for gene specific information),
#' \code{\link{assays}} (for gene by cell matrices), 
#' or \code{\link{metadata}} (for parameters configurations) slots. 
#'
#'
#' @return SingleCellExperiment object containing the simulated counts and
#' intermediate values.
#'
#' @references
#' Tian J, Wang J, Roeder K. ESCO: single cell expression simulation incorporating gene co-expression. 
#' bioRxiv. 2020.
#'
#' Paper: \url{https://www.biorxiv.org/content/10.1101/2020.10.20.347211v1}
#'
#' Code: \url{https://github.com/JINJINT/}
#'
#' @seealso
#' \code{\link{escoSimLibSizes}}, 
#' \code{\link{escoSimGeneMeans}},
#' \code{\link{escoSimMeans}},
#' \code{\link{escoSimDE}}, 
#' \code{\link{escoSimCellMeans}}, 
#' \code{\link{escoSimTrueCounts}},
#' \code{\link{escoSimDropout}}
#'
#' @examples
#' # Simulation with default parameters
#' sim <- escoSimulate()
#' \dontrun{
#' # Simulation with different number of genes
#' sim <- escoSimulate(nGenes = 1000)
#' # Simulation with custom parameters
#' params <- newescoParams(nGenes = 100, mean.rate = 0.5)
#' sim <- escoSimulate(params)
#' # Simulation with adjusted custom parameters
#' sim <- escoSimulate(params, mean.rate = 0.6, out.prob = 0.2)
#' # Simulate group
#' sim <- escoSimulate(type = "group")
#' # Simulate tree
#' sim <- escoSimulate(type = "tree")
#' # Simulate traj
#' sim <- escoSimulate(type = "traj")
#' }
#' @importFrom SummarizedExperiment rowData colData colData<- assays 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods validObject
#' @importFrom svMisc progress
#' @export
#' @rdname escoSimulate

escoSimulate <- function(params = newescoParams(), 
                          type = c("single", "group", "tree", "traj"), 
                          verbose = TRUE, numCores=2, ...){

    checkmate::assertClass(params, "escoParams")

    type <- match.arg(type)

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
    nGroups <- getParam(params, "nGroups")
    group.prob <- getParam(params, "group.prob")
    dirname <- getParam(params, "dirname")
    trials <- getParam(params, "trials")
    
    if(type == "tree"){
      tree <- getParam(params, "tree")[[1]]
      nG <- length(tree$tip.label)
      if(nG!=nGroups){
        message("The given number of Groups and Tree structure does not match, 
                thus recalibarting to tree strucutrue with equal size groups....")
        group.prob = rep(round(1/nG,3), nG)
        group.prob[1] = 1 - sum(group.prob[-1])
        params <- setParams(params, group.prob = group.prob)
      }
    }

    if (nGroups == 1 && type %in% c("group", "tree")) {
        warning("nGroups is 1, switching to single mode")
        type <- "single"
    }

    if (verbose) {message("Creating simulation object...")}
    
    # Set up name vectors
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    if (type %in% c("group", "tree")){
        group.names <- paste0("Group", seq_len(nGroups))
    }
    
    # Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    if (type %in% c("group", "tree")) {
        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- group.names[groups]
    }

    if (verbose) {message("Simulating library sizes...")}
    sim <- escoSimLibSizes(sim, verbose)
    
    if (verbose) {message("Simulating gene means...")}
    sim <- escoSimGeneMeans(sim, verbose)
    
    if (type == "single") {
        sim <- escoSimSingleCellMeans(sim, verbose)
    }
    if (type == "group") {
        if (verbose) {message("Simulating group DE...")}
        sim <- escoSimGroupDE(sim, verbose)
        if (verbose) {message("Simulating cell means...")}
        sim <- escoSimGroupCellMeans(sim, verbose)
    } 
    if(type == "tree") { 
      if (verbose) {message("Simulating tree structure...")}
      sim <- escoSimTreeDE(sim, verbose)
      if (verbose) {message("Simulating cell means...")}
      sim <- escoSimTreeCellMeans(sim, verbose)
    }
    if(type == "traj") { 
      if (verbose) {message("Simulating trajectory...")}
      sim <- escoSimPathDE(sim, verbose)
      if (verbose) {message("Simulating cell means...")}
      sim <- escoSimPathCellMeans(sim, verbose)
    }
    
    if (verbose) {message("Simulating true counts...")}
    sim <- escoSimTrueCounts(sim, type, verbose,numCores)
    
    
    params <- metadata(sim)$Params
    trials <- getParam(params, "trials")

    if(trials>1 & (!dir.exists(dirname))){
        message("Warning: detect calls for multiple trials, but no directory to save files...")
        message("Modify: only do the simulations with one trial...")
        params <- setParams(params, trials = 1)
        metadata(sim)$Params= params
    }
    
    dropout.type = getParam(params, "dropout.type")
    
    if(!dir.exists(dirname)){
      if (verbose) message("Adding technical noise ...")
      if(length(dropout.type)>0){
        if("zeroinflate" %in% dropout.type)sim<-Observed_Counts(sim, "", verbose,numCores) 
        if("downsample" %in% dropout.type)sim<-escoSimDropout(sim, "", verbose)
      }
      return(sim)
    }
    else{
      if (verbose) message("Adding technical noise ...")
      trials<-getParam(metadata(sim)$Params, "trials")
      
      if(trials>1){
        if(verbose)message("Starting multiple trials.....")
        for(trial in 1:trials){
          if(dropout.type!="none"){
            if("downsample" %in% dropout.type)sim<-Observed_Counts(sim, trial,verbose,numCores) 
            if("zeroinflate" %in% dropout.type)sim<-escoSimDropout(sim, trial,verbose)
          }
        }
      }
      else{
        if(dropout.type!="none"){
          if("downsample" %in% dropout.type)sim<-Observed_Counts(sim, "",verbose,numCores) 
          if("zeroinflate" %in% dropout.type)sim<-escoSimDropout(sim, "",verbose)
        }
      }
      if(dir.exists(dirname))saveRDS(sim, paste0(dirname, "sim.rds"))
      return(sim)
    }
    

}

#' @rdname escoSimulate
#' @export
escoSimulateSingle <- function(params = newescoParams(), 
                                verbose = TRUE, ...) {
    sim <- escoSimulate(params = params, type = "single", 
                         verbose = verbose, ...)
    return(sim)
}

#' @rdname escoSimulate
#' @export
escoSimulateGroups <- function(params = newescoParams(), 
                                verbose = TRUE, ...) {
    sim <- escoSimulate(params = params, type = "group", 
                         verbose = verbose, ...)
    return(sim)
}

#' @rdname escoSimulate
#' @export
escoSimulateTree <- function(params = newescoParams(), 
                                verbose = TRUE, ...) {
  sim <- escoSimulate(params = params, type = "tree", 
                       verbose = verbose, ...)
  return(sim)
}

#' @rdname escoSimulate
#' @export
escoSimulateTraj <- function(params = newescoParams(), 
                              verbose = TRUE, ...) {
  sim <- escoSimulate(params = params, type = "traj", 
                      verbose = verbose, ...)
  return(sim)
}


#' Simulate library sizes
#'
#' Simulate expected library sizes. Typically a log-normal distribution is used
#' but there is also the option to use a normal distribution. In this case any
#' negative values are set to half the minimum non-zero value.
#'
#' @param sim SingleCellExperiment to add library size to.
#'
#' @return SingleCellExperiment with simulated library sizes.
#'
#' @importFrom SummarizedExperiment colData colData<- 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats rlnorm rnorm
#' @rdname escoSimLib
escoSimLibSizes <- function(sim, verbose) {
    params <- metadata(sim)$Params
    nCells <- getParam(params, "nCells")
    Groups <- colData(sim)$Group
    lib.loc <- getParam(params, "lib.loc")
    lib.scale <- getParam(params, "lib.scale")
    lib.norm <- getParam(params, "lib.norm")
    lib.method <- getParam(params, "lib.method")
    
    if(lib.method == "fit"){
      if(verbose)message("Sampling from fitted log normal / normal...")
      if(length(lib.loc)>1){
        colData(sim)$ExpLibSize  <- rep(0,nCells)
        for(idx in 1:length(lib.loc)){
          cells = which(colData(sim)$Group==paste0("Group",idx))
          if (lib.norm[idx]) {
            exp.lib.sizes <- rnorm(length(cells), lib.loc[idx], lib.scale[idx])
            min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
            exp.lib.sizes[exp.lib.sizes < 0] <- min.lib / 2
          } else {
            exp.lib.sizes <- rlnorm(cells, lib.loc[idx], lib.scale[idx])
          }
          colData(sim)$ExpLibSize[cells] = exp.lib.sizes
        }
      }
      else{
        if (lib.norm) {
          exp.lib.sizes <- rnorm(nCells, lib.loc, lib.scale)
          min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
          exp.lib.sizes[exp.lib.sizes < 0] <- min.lib / 2
        } else {
          exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
        }
      }
    }
    
    if(lib.method == "density"){
      if(verbose)message("Sampling from provided density object...")
      lib.dens <- getParam(params, "lib.dens")
      exp.lib.sizes <- sampleDensity(nCells, lib.dens)
    }
    colData(sim)$ExpLibSize  <- exp.lib.sizes
    
    metadata(sim)$Params = params
    
    return(sim)
}

#' Simulate gene means
#'
#' Simulate gene means from a gamma distribution. Also simulates outlier
#' expression factors. Genes with an outlier factor not equal to 1 are replaced
#' with the median mean expression multiplied by the outlier factor.
#'
#' @param sim SingleCellExperiment to add gene means to.
#'
#' @return SingleCellExperiment with simulated gene means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats rgamma median
#' @rdname escoSimMean
escoSimGeneMeans <- function(sim, verbose) {
    params <- metadata(sim)$Params
    
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")
    mean.method <- getParam(params, "mean.method")
    mean.dens <- getParam(params, "mean.dens")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene

    # Simulate base gene means
    if(mean.method =="density"){
      if(verbose)message("Sampling from provided density object...")
      mean.dens <- getParam(params, "mean.dens")
      base.means.gene <- exp(sampleDensity(nGenes, mean.dens, lower = -Inf))
      
      }
    if(mean.method =="fit"){
      if(verbose)message("Sampling from fitted gamma...")
      base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
      
    }

    # Add expression outliers
    is.selected = as.logical(rbinom(nGenes, 1, out.prob))
    outlier.facs <- getLNormFactors(nGenes, is.selected, 0, out.facLoc, out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs
    rowData(sim)$GeneMean <- means.gene
    
    batch.means.cell <- matrix(1, ncol = nCells, nrow = nGenes) * means.gene
    colnames(batch.means.cell) <- cell.names
    rownames(batch.means.cell) <- gene.names
    assays(sim)$Means <- batch.means.cell
    
    metadata(sim)$Params = params
    
    return(sim)
}

#' Simulate group differential expression
#'
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taked to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#'
#' @return SingleCellExperiment with simulated differential expression.
#'
#' @rdname escoSimDE
#' @importFrom SummarizedExperiment rowData 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats ecdf
escoSimGroupDE <- function(sim, verbose) {
    params <- metadata(sim)$Params
    
    nGenes <- getParam(params, "nGenes")
    
    gene.names <- rowData(sim)$Gene
    gene.means <- rowData(sim)$GeneMean
    #percen <- ecdf(gene.means)
    #percentile <- percen(gene.means)
    
    de.prob <- getParam(params, "de.prob")
    nGroups <- getParam(params, "nGroups")
    #de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")
    deall.prob <- getParam(params, "deall.prob")
    #de.quantile <-  getParam(params, "de.quantile")
    corr<-getParam(params, "corr")
    de.rank <- getParam(params, "de.rank")
    de.facrank <- getParam(params, "de.facrank")
    
    if(length(corr)==nGroups+1){
      len = lapply(corr[2:length(corr)], nrow)
      len = unlist(len, use.names=FALSE)
      # for(idx in 1:length(de.quantile)){
      #   which(percentile == de.quantile[idx])
      #   deall.genes = sample(percentile, len[idx])
      # }
      deall.genes = sample(nGenes, sum(len))
      deall.prob = length(deall.genes)/nGenes
      de.prob = len/length(deall.genes)
      params<-setParams(params, deall.prob = deall.prob, de.prob = de.prob)
      marker.genes = rep(1:nGroups, len)
      markers = rep(0, nGenes)
      markers[deall.genes] = marker.genes
    }
    else{
      if(length(de.rank)==nGroups){
        markers = rep(0, nGenes)
        ordsim = sort(gene.means, index.return = TRUE)$ix
        for(idx in 1:nGroups){
          markers[which(ordsim%in%de.rank[[idx]])] = idx
        }
        deall.genes = which(markers!=0)
        rowData(sim)$DEgenes = rep(0, nGenes)
        rowData(sim)$DEgenes[deall.genes] = 1
        deall.prob = length(deall.genes)/nGenes
        len = lapply(de.rank, length)
        len = unlist(len, use.names=FALSE)
        de.prob = len/length(deall.genes)
        params<-setParams(params, deall.prob = deall.prob, de.prob = de.prob)
      }
      else{
        deall.genes = sample(nGenes, floor(deall.prob*nGenes))
        rowData(sim)$DEgenes = rep(0, nGenes)
        rowData(sim)$DEgenes[deall.genes] = 1
        marker.genes = sample(seq_len(nGroups), size = length(deall.genes), replace = TRUE, prob = de.prob)
        markers = rep(0, nGenes)
        markers[deall.genes] = marker.genes
      }
    }
      
    for(idx in seq_len(nGroups)){
        markershere = (markers==idx)
        de.facs <- getLNormFactors(nGenes, markershere, 0,
                                   de.facLoc[idx], de.facScale[idx])
        if(length(de.rank)==nGroups){
          ordsim = sort(gene.means, index.return = TRUE)$ix
          ord = match(ordsim[markershere], de.rank[[idx]])
          ordfac = de.facrank[[idx]][ord]
          de.facs[markershere] = sort(de.facs[markershere])[ordfac]
        }
        rowData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
    }
    
    rowData(sim)$GeneGroup = markers
    
    metadata(sim)$Params = params

    return(sim)
}


#' @rdname escoSimDE
#' @importFrom SummarizedExperiment rowData 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom ape vcv.phylo
#' @importFrom MASS mvrnorm
escoSimTreeDE <- function(sim, verbose) {
  params <- metadata(sim)$Params
  
  tree<- getParam(params, "tree")[[1]]
  de.center<-getParam(params, "de.center")
  deall.prob <- getParam(params, "deall.prob")
  de.prob <- getParam(params, "de.prob")
  
  nGenes <- getParam(params, "nGenes")
  nGroups <- getParam(params, "nGroups")
  corr<-getParam(params, "corr")
  
  if(length(corr)>0 & length(corr)!=nGroups){
    params<-setParams(params, corr = list())
    cat("Sorry, cannot use given correlation matrix when simulating tree strucutre.....\n
        Using randmized instead...\n")
  }
  
  gene.names <- rowData(sim)$Gene
  
  de.prob <- getParam(params, "de.prob")
  de.downProb <- getParam(params, "de.downProb")
  de.facLoc <- getParam(params, "de.facLoc")
  de.facScale <- getParam(params, "de.facScale")
  
  deall.genes = sample(nGenes, floor(deall.prob*nGenes))
  rowData(sim)$DEgenes = rep(0, nGenes)
  rowData(sim)$DEgenes[deall.genes] = 1
  marker.genes = sample(seq_len(nGroups), size = length(deall.genes), replace = TRUE, prob = de.prob)
  markers = rep(0, nGenes)

  markers[deall.genes] = 1
  treecor <- vcv.phylo(tree, cor=T)
  treecor = sqrt(treecor)
  de.facLocvec = matrix(0, nGenes, nGroups)
  de.facLocvec[deall.genes,] <- mvrnorm(length(deall.genes), rep(de.center, nGroups), treecor)
  for(idx in seq_len(nGroups)){
      markers = (markers>0)
      de.facs <- getLNormFactors(nGenes, markers, 0,
                                 de.facLocvec[,idx], de.facScale[idx])
      rowData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
  }
  metadata(sim)$Params = params
  return(sim)
}


#' Simulate trajectory
#'
#' Simulate gene means for each step along each path of a esco simulation
#'
#' @param params escoParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @importFrom S4Vectors metadata metadata<-
#'
#' @details
#' The method of simulating paths is inspired by the method used in the PROSSTT
#' simulation. Changes in expression are controlled by \code{paths.nPrograms}
#' regulatory programs. Each of the regulatory genes in the gene network has
#' some association with each program. This is analogous to there being changes
#' in the environment (the programs) which are sensed by receptors (regulatory
#' genes) and cause changes in expression downstream. For each path a random
#' walk is generated for each program and the changes passed on to the
#' regulatory genes. At each step the changes propagate through the network
#' according to the weights on edges between genes. This algorithm is fairly
#' simple but should result in correlation relationships between genes. However
#' it is likely to be improved and adjusted in the future.
#'
#' The path structure itself is specified by the \code{paths.design} parameter.
#' This is a \code{data.frame} with three columns: "traj", "From", and "Steps".
#' The Path field is an ID for each path while the Steps field controls the
#' length of each path. Increasing the number of steps will increase the
#' difference in expression between the ends of the paths. The From field sets
#' the originating point of each path. For example a From of \code{0, 0, 0}
#' would indicate three paths from the origin while a From of \code{0, 1, 1}
#' would give a branching structure with Path 1 beginning at the origin and
#' Path 2 and Path 3 beginning at the end of Path 1.
#'
#' @references
#'
#' Papadopoulos N, Parra RG, Söding J. PROSSTT: probabilistic simulation of
#' single-cell RNA-seq data for complex differentiation processes.
#' Bioinformatics (2019). \url{https://doi.org/10.1093/bioinformatics/btz078}.
#'
#' @return SingleCellExperiment with simulated trajetory means.
#' @rdname escoSimDE
escoSimPathDE <- function(sim, verbose) {
  params = metadata(sim)$Params
  if (verbose) {message("Simulating trajetories...")}
  nGenes <- getParam(params, "nGenes")
  corr <- getParam(params, "corr")
  withcorr <- getParam(params, "withcorr")
  
  paths.design <- getParam(params, "paths.design")

  if(withcorr & (length(corr)>0)){
    deall.genes = sample(nGenes, nrow(corr[[1]]))}
  else{
    paths.deprob <- getParam(params, "paths.deprob")
    deall.genes = sample(nGenes, floor(paths.deprob*nGenes))
  }
  paths.DEgenes = rep(0, nGenes)
  paths.DEgenes[deall.genes] = 1
  
  paths.changes <- vector("list", nrow(paths.design))
  paths.velocity <- vector("list", nrow(paths.design))
  paths.factors <- vector("list", nrow(paths.design))
  paths.diffusion <- vector("list", nrow(paths.design))
  
  paths.graph <- igraph::graph_from_data_frame(paths.design)
  paths.order <- names(igraph::topo_sort(paths.graph, mode = "in"))
  paths.order <- as.numeric(paths.order)
  
  # Remove the origin because it is not a path
  paths.order <- paths.order[paths.order != 0]
  
  for (path in paths.order) {
    if (verbose) {message("Simulating path ", path, "...")}
    nSteps <- paths.design$Steps[path]
    from <- paths.design$From[path]
    
    changes <- matrix(0, nrow = nGenes, ncol = nSteps + 1)
    velocity <- matrix(0, nrow = length(deall.genes), ncol = nSteps + 1)
    diffusion <- matrix(0, nrow = length(deall.genes), ncol = nSteps + 1)
    if(from==0){
      velocity[,1] = rnorm(length(deall.genes), mean = 0, sd = 0.02)
      diffusion[,1] = rep(0, length(deall.genes))
    }
    else{
      velocity[,1] = rnorm(length(deall.genes), mean = 0, sd = 0.02)
      diffusion[,1] = paths.diffusion[[from]][, ncol(paths.diffusion[[from]])]
    }
    
    
    for (step in seq_len(nSteps) + 1) {
      
      velocity[,step] <-  velocity[, step - 1] +  rnorm(deall.genes, mean = 0, sd = 2/(nSteps))
      diffusion[,step] <- diffusion[,step-1] + velocity[, step - 1]
      velocity[,step] <- sapply(velocity[, step], function(x)sign(x)*min(abs(x), 6/(nSteps)))
      changes[deall.genes, step] <- diffusion[, step]
    }
    
    if (from == 0) {
      factors <- changes[, seq_len(nSteps)]
      
    } else {
      factors <- changes[, seq_len(nSteps) + 1]
    }
    paths.changes[[path]] <- changes
    paths.factors[[path]] <- factors
    paths.diffusion[[path]] <- diffusion
    paths.velocity[[path]] <- velocity
  }
  
  params <- setParams(params, paths.factors = paths.factors, paths.DEgenes = paths.DEgenes)
  metadata(sim)$Params<-params
  return(sim)
}



#' Simulate cell means
#'
#' Simulate a gene by cell matrix giving the mean expression for each gene in
#' each cell. Cells start with the mean expression for the group they belong to
#' (when simulating groups) or cells are assigned the mean expression from a
#' random position on the appropriate path (when simulating tree). The selected
#' means are adjusted for each cell's expected library size.
#'
#' @param sim SingleCellExperiment to add cell means to.
#'
#' @return SingleCellExperiment with added cell means.
#'
#' @rdname escoSimCellMeans
#' @importFrom SummarizedExperiment rowData colData assays assays<- 
#' @importFrom S4Vectors metadata metadata<-
escoSimSingleCellMeans <- function(sim, verbose) {
    params <- metadata(sim)$Params
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$Means

    cell.means.gene <- batch.means.cell
    cell.props.gene <- t(t(cell.means.gene)/colSums(cell.means.gene))

    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
    
    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell
    
    metadata(sim)$Params = params
    return(sim)
}

#' @importFrom SummarizedExperiment rowData colData assays assays<- 
#' @importFrom S4Vectors metadata metadata<-
#' @rdname escoSimCellMeans
escoSimGroupCellMeans <- function(sim, verbose) {
    params <- metadata(sim)$Params
    nGenes <- getParam(params, "nGenes")
    de.prob <- getParam(params, "de.prob")
    de.facLoc <- getParam(params, "de.facLoc")
    nGroups <- getParam(params, "nGroups")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    groups <- colData(sim)$Group
    group.names <- sort(unique(groups))
    lib.loc <- getParam(params, "lib.loc")
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$Means
    

    group.facs.gene <- rowData(sim)[, paste0("DEFac", group.names)]
    DEgene.name = gene.names[which(rowData(sim)$DEgenes==1)]
    cell.facs.gene <- as.matrix(group.facs.gene[, paste0("DEFac", groups)])
    cell.means.gene <- batch.means.cell*cell.facs.gene
  
    cell.props.gene <- t(t(cell.means.gene)/colSums(cell.means.gene))
    
    base.means.cell <- t(t(cell.props.gene)*exp.lib.sizes)
    

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell
    
    metadata(sim)$Params = params
    return(sim)
}

#' @rdname escoSimCellMeans
#' @importFrom SummarizedExperiment rowData colData assays assays<- 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SC3 get_marker_genes
#' @importFrom doBy which.minn which.maxn
escoSimTreeCellMeans <- function(sim, verbose) {
  params<-metadata(sim)$Params
  nGenes <- getParam(params, "nGenes")
  de.prob <- getParam(params, "de.prob")
  de.facLoc <- getParam(params, "de.facLoc")
  nGroups <- getParam(params, "nGroups")
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  groups <- colData(sim)$Group
  group.names <- sort(unique(groups))
  exp.lib.sizes <- colData(sim)$ExpLibSize
  batch.means.cell <- assays(sim)$Means
  group.facs.gene <- rowData(sim)[, paste0("DEFac", group.names)]
  DEgene.name = gene.names[which(rowData(sim)$DEgenes==1)]
  cell.facs.gene <- as.matrix(group.facs.gene[, paste0("DEFac", groups)])
  cell.means.gene <- batch.means.cell*cell.facs.gene
  
  # identify the marker genes
  markers = get_marker_genes(cell.means.gene[as.character(DEgene.name),],
                               as.factor(groups))
  GeneGroup = rep(0, nGenes)
  mulmark = matrix(1, nGenes, ncol(cell.means.gene))
  for(i in 1:length(as.vector(table(markers$clusts)))){
      index = which(markers$clusts == i)
      pval = markers$pval[index]
      len = as.integer(de.prob[i]*length(index))
      GeneGroup[index[which.minn(pval, len)]] = i
      mulmark[index[which.minn(pval, len)], which(groups==paste0('Group',i))] = de.facLoc[i]
  }
  rowData(sim)$GeneGroup = GeneGroup
  cell.means.gene = cell.means.gene * mulmark
  
  cell.props.gene <- t(t(cell.means.gene)/colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene)*exp.lib.sizes)
  
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  assays(sim)$BaseCellMeans <- base.means.cell
  
  metadata(sim)$Params = params
  return(sim)
}



#' Simulate esco trajetory cell means
#'
#' @param sim SingleCellExperiment containing simulation.
#' @param params escoParams object with simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @importFrom S4Vectors metadata metadata<-
#'
#' @details
#' Cells are first assigned to a path and a step along that path. This is
#' controlled by the \code{cells.design} parameter which is a \code{data.frame}
#' with the columns "traj", "Probability", "Alpha" and "Beta". The Path field
#' is an ID for each path and the Probability field is the probability that a
#' cell will come from that path (must sum to 1). The Alpha and Beta parameters
#' control the density of cells along the path. After they are assigned to paths
#' the step for each cell is sampled from a Beta distribution with parameters
#' shape1 equals Alpha and shape2 equals beta. The distribution
#' can be viewed using \code{hist(rbeta(10000, Alpha, Beta), breaks = 100)}.
#' Once cells are assigned to paths and steps the correct means are extracted
#' from the \code{paths.means} parameter and adjusted based on each cell's
#' library size. An adjustment for BCV is then applied. 
#'
#' @return SingleCellExperiment with cell means
#' @rdname escoSimCellMeans
escoSimPathCellMeans <- function(sim, verbose) {
  params = metadata(sim)$Params
  
  cell.names <- colData(sim)$Cell
  gene.names <- rowData(sim)$Gene
  nGenes <- getParam(params, "nGenes")
  nCells <- getParam(params, "nCells")
  cells.design <- getParam(params, "cells.design")
  paths.design <- getParam(params, "paths.design")
  paths.factors <- getParam(params, "paths.factors")
  genemeans <- rowData(sim)$GeneMean
  
  paths.means <- lapply(paths.factors, function(x) {
    ans = 2^x * genemeans
    return(ans)
  })
  names(paths.means) <- paste0("traj", paths.design$Path)
  params <- setParams(params, paths.means = paths.means)
  
  exp.lib.sizes <- colData(sim)$ExpLibSize
  
  if (verbose) {message("Assigning cells to paths...")}
  cells.paths <- sample(cells.design$Path, nCells, replace = TRUE,
                        prob = cells.design$Probability)
  
  if (verbose) {message("Assigning cells to steps...")}
  paths.cells.design <- merge(paths.design, cells.design)
  steps.probs <- apply(paths.cells.design, 1, function(path) {
    steps <- path["Steps"]
    probs <- getBetaStepProbs(path["Steps"], path["Alpha"], path["Beta"])
    
    # Return a list to avoid getting a matrix if all path lengths are equal
    return(list(probs))
  })
  # Remove unnecessary list level
  steps.probs <- lapply(steps.probs, "[[", 1)
  names(steps.probs) <- paths.cells.design$Path
  
  cells.steps <- vapply(cells.paths, function(path) {
    probs <- steps.probs[[path]]
    step <- sample(seq_len(length(probs)), 1, prob = probs)
    return(step)
  }, c(Step = 0))
  
  if (verbose) {message("Simulating cell means...")}
  basecells.means <- vapply(seq_len(nCells), function(cell) {
    path <- cells.paths[cell]
    step <- cells.steps[cell]
    means <- paths.means[[path]][, step]
    return(means)
  }, as.numeric(seq_len(nGenes)))

  # Adjust mean based on library size
  basecells.means <- exp.lib.sizes * t(t(basecells.means) / colSums(basecells.means))

  colnames(basecells.means) <- cell.names
  rownames(basecells.means) <- gene.names
  assays(sim)$BaseCellMeans <- basecells.means
  
  colData(sim)$Path = cells.paths
  colData(sim)$Step = cells.steps  
  metadata(sim)$Params<-params
  return(sim)
}



#' Simulate true counts
#'
#' Simulate a true counts matrix. Counts are simulated from a poisson
#' distribution where Each gene in each cell has it's own mean based on the
#' group (or path position), expected library size and BCV.
#'
#' @param sim SingleCellExperiment to add true counts to.
#'
#' @return SingleCellExperiment with simulated true counts.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<- 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats rpois qnbinom pnorm cor
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import progress
#' @rdname escoSimTruth
escoSimTrueCounts <- function(sim, type, verbose, numCores = 2) {
    params<-metadata(sim)$Params
    withcorr <- getParam(params, "withcorr")
    corr<-getParam(params, "corr")
    
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    bcv.common <- getParam(params, "bcv.common")
    bcv.df <- getParam(params, "bcv.df")
    basecell.means <- assays(sim)$BaseCellMeans
    basecell.means <- as.matrix(basecell.means)
    dirname <- getParam(params, "dirname")

    house.prob <- getParam(params, "house.prob")
    nGroups <- getParam(params, "nGroups")
    
    
    if(type %in% c("group", "tree")){
      groups = colData(sim)$Group
    }
    if(type=="single"){
      groups = rep(1, nCells)
    }

    bcv <- (bcv.common + (1 / sqrt(basecell.means)))*sqrt(bcv.df/rchisq(nGenes, df = bcv.df))
    dimnames(bcv) =  dimnames(basecell.means)
    bcv = as.matrix(bcv)
    
    cell.means <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                                scale = basecell.means * (bcv ^ 2)),
                         nrow = nGenes, ncol = nCells)
    
    # use gaussian copular to simulate co-expression
    if(!withcorr){
      true.count <- matrix(rpois(as.numeric(nGenes) * as.numeric(nCells),
                                 lambda = cell.means),
                           nrow = nGenes, ncol = nCells)
    }
    
    else{ 
      
      if(type == "single"){
        if(length(corr)!=0){
          corrgenes = sort(sample(seq_len(nGenes), nrow(corr[[1]])))
          genemeans = rowData(sim)$GeneMean
          ordsim = sort(genemeans[corrgenes], index.return = TRUE)$ix
          rho = makespd(corr[[1]][ordsim, ordsim])
          rownames(rho) = gene.names[corrgenes]
          colnames(rho) = rownames(rho)  
        }
        else{
          corr.prob = getParams(params, "corr.prob")$corr.prob
          corrgenes = sample(seq_len(nGenes), as.integer(corr.prob*nGenes))
          rho = randcor(length(corrgenes))
          rownames(rho) = gene.names[corrgenes]
          colnames(rho) = rownames(rho)  
        }
        params<-setParams(params, corr = list(rho))
        copular = randcop(rho, nCells)
        
        mulpo<-function(i){
          count = rnbinom(nGenes, size = 1/(bcv[,i]^2), mu = basecell.means[,i])
          count[corrgenes] = qnbinom(copular[,i],  size = 1/(bcv[corrgenes,i]^2), mu = basecell.means[corrgenes,i])
          return(count)
        }
      } 
      
      if(type %in% c("group","tree")){
        housegenes = gene.names
        for(id in sort(unique(groups))){
          DEindex = which(rowData(sim)[[paste0("DEFac", id)]]!=1)
          housegenes[DEindex] = NA 
        }
        housegenes = housegenes[!is.na(housegenes)]
        housegenes = as.character(housegenes)
        if(length(corr)==nGroups+1)housegenes = housegenes[which.maxn(rowMeans(basecell.means)[housegenes], nrow(corr[[1]]))]
        else housegenes = housegenes[which.maxn(rowMeans(basecell.means)[housegenes], ceiling(house.prob*length(housegenes)))]
        housegenes = as.character(housegenes)
        
        if(length(corr)!=0){
          if(length(corr)!=nGroups+1){
            cat("The given list of correlation matrix should be of length 1 + number of groups...\n
                 Using the randomized instead...\n")
          }
          else{
            rholist = list()
            genemeans = rowData(sim)$GeneMean
            ordsim = sort(genemeans[match(housegenes, gene.names)], index.return  = TRUE)$ix
            rholist[["housekeep"]] = makespd(corr[[1]][ordsim, ordsim])
            rownames(rholist[["housekeep"]]) = housegenes
            colnames(rholist[["housekeep"]]) = rownames(rholist[["housekeep"]])
            for(idx in 1:length(unique(groups))){
                ordsim = sort(genemeans[which(rowData(sim)$GeneGroup==idx)], index.return  = TRUE)$ix
                rholist[[paste0("Group",idx)]] = makespd(corr[[idx+1]][ordsim, ordsim])
                rownames(rholist[[paste0("Group",idx)]]) = gene.names[which(rowData(sim)$GeneGroup==idx)]
                colnames(rholist[[paste0("Group",idx)]]) = rownames(rholist[[paste0("Group",idx)]])
            }
          }
        }
        else{
          rholist = list()
          rholist[["housekeep"]] = randcor(length(housegenes))
          rownames(rholist[["housekeep"]]) = housegenes
          colnames(rholist[["housekeep"]]) = rownames(rholist[["housekeep"]])
          for(idx in 1:length((unique(groups)))){
              markgenes = as.character(gene.names[which(rowData(sim)$GeneGroup==idx)])
              rholist[[paste0("Group",idx)]] = randcor(length(markgenes))
              rownames(rholist[[paste0("Group",idx)]]) = markgenes
              colnames(rholist[[paste0("Group",idx)]]) = rownames(rholist[[paste0("Group",idx)]])
          }
        }
        params<-setParams(params, corr = rholist)
        coplist = lapply(rholist, function(covmat)randcop(covmat, nCells))
        
        mulpo<-function(i){
          count = rnbinom(nGenes, size = 1/(bcv[,i]^2), mu = basecell.means[,i])
          names(count) = gene.names 
          
          idx = which(sort(unique(groups))==groups[i])
          DEgenes = as.character(gene.names[which(rowData(sim)$GeneGroup==idx)])
          copular = coplist[[idx+1]]
          count[DEgenes] = qnbinom(copular[,i], size = 1/(bcv[DEgenes,i]^2), mu = basecell.means[DEgenes,i])
          
          housecopular = coplist[[1]]
          count[housegenes] = qnbinom(housecopular[,i], size = 1/(bcv[housegenes,i]^2), mu = basecell.means[housegenes,i])
          
          return(count)
        }
      }
      
      if(type=="traj"){
        paths.DEgenes = getParam(params, "paths.DEgenes")
        corrgenes = which(paths.DEgenes==1)
        
        if(length(corr)>0){
          rho = corr[[1]]
        }
        else{
          rho = randcor(length(corrgenes))
        }
        
        rownames(rho) = gene.names[corrgenes]
        colnames(rho) = rownames(rho)  
        
        params<-setParams(params, corr = list(rho))
        copular = randcop(rho, nCells)
        
        mulpo<-function(i){
          count = rnbinom(nGenes, size = 1/(bcv[,i]^2), mu = basecell.means[,i])
          count[corrgenes] = qnbinom(copular[,i],  size = 1/(bcv[corrgenes,i]^2), mu = basecell.means[corrgenes,i])
          return(count)
        }
        
      }
      
      # parallel
    
      if (is.null(numCores)){
        numCores <- parallel::detectCores()
      }
      cl <- makeCluster(numCores)
      registerDoSNOW(cl)
      total <- ncol(basecell.means)
      
      pb <- progress_bar$new(
        format = "progress = :letter [:bar] :elapsed | eta: :eta", total = total, width = 60)
      
      progress <- function(n){
        pb$tick(tokens = list(letter = rep("", total)[n]))
      } 
      opts <- list(progress = progress)
      true.count <- foreach(i = 1:total, .combine = cbind, .options.snow = opts, .export = c("rowData")) %dopar% {
        return(mulpo(i))
      }
      stopCluster(cl) 
    }
    
    colnames(true.count) <- cell.names
    rownames(true.count) <- gene.names
    if(dir.exists(dirname)){
      if(verbose)cat("Saving true counts: ", 100*round(mean(true.count == 0), 4), "% zeros...\n")
      saveRDS(true.count, paste0(dirname,"True.rds"))
    }
    if(sum(true.count==Inf)>0 & verbose)cat("Detecting",100*round(sum(true.count==Inf)/(nCells*nGenes),4),"% Inf true count value, modifying it to max + 10....\n")
    true.count[true.count==Inf] = max(true.count[true.count!=Inf]) + 10
    assays(sim)$TrueCounts <- true.count
    assays(sim)$CellMeans <- cell.means
    
    metadata(sim)$Params = params
    return(sim)
}

#' Simulate dropout
#'
#' A logistic function is used to form a relationshop between the expression
#' level of a gene and the probability of dropout, giving a probability for each
#' gene in each cell. These probabilities are used in a Bernoulli distribution
#' to decide which counts should be dropped.
#' @rdname escoSimNoise
#' @param sim SingleCellExperiment to add dropout to.
#'
#' @return SingleCellExperiment with simulated dropout and observed counts.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats rbinom
#' @export
escoSimDropout <- function(sim, trial, verbose) {
    params <- metadata(sim)$Params
    dirname <- getParam(params, "dirname")
    true.counts <- assays(sim)$TrueCounts
   
    dropout.mid <- getParam(params, "dropout.mid")
    dropout.shape <- getParam(params, "dropout.shape")
    dropout.cort <- getParam(params, "dropout.cort")
    
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    bcv.df <- getParam(params, "bcv.df")
    bcv.common <- getParam(params, "bcv.common")
    nGenes <- getParam(params, "nGenes")
    nGroups <- getParam(params, "nGroups")
    withcorr <- getParam(params, "withcorr")
    genemeans <- rowData(sim)$GeneMean
    
    if(withcorr)cell.means <- assays(sim)$CellMeans
    else cell.means <- assays(sim)$BaseCellMeans
    
    cell.normmeans  = median(colSums(cell.means)) * t(t(cell.means)/colSums(cell.means))
      
    if(length(dropout.mid)>1 & (!dir.exists(dirname))){
      message("Warning: detect calls for simulating multiple configuration, but no directory to save files....")
      message("Modifying: only do simulation for the first configuration....")
      dropout.mid = dropout.mid[1]
    }
    
    for(dropout.midd in dropout.mid){
      dropout.shape <- getParam(params, "dropout.shape")
      dm = dropout.midd
      dropout.midd <- rep(dropout.midd, nCells)
      dropout.shape <- rep(dropout.shape, nCells)
      
      # Generate probabilites based on expression
      drop.prob <- sapply(seq_len(nCells), function(idx) {
          eta <- log(cell.normmeans[,idx])
          return(logistic(eta, x0 = dropout.midd[idx], k = dropout.shape[idx]))
      })
        
      if(!dropout.cort)keep.prob <- 1 - drop.prob 
      else{
          keep.prob  <- (1 - drop.prob)/(1-dpois(0, lambda = cell.normmeans))
          # bcv <- (bcv.common + (1 / sqrt(cell.means)))*sqrt(bcv.df/rchisq(nGenes, df = bcv.df))
          # dimnames(bcv) =  dimnames(cell.means)
          # bcv = as.matrix(bcv)
          # keep.prob <- (1-drop.prob)/(1-dnbinom(0, size = 1/(bcv^2), mu = cell.means))
        }
        
      keep.prob[keep.prob>1] <- 1
      keep.prob[keep.prob<0] <- 0
      keep.prob[is.na(keep.prob)] <- 1
        
      keep <- matrix(rbinom(nCells * nGenes, 1, keep.prob),
                       nrow = nGenes, ncol = nCells)

        
      counts <- true.counts * keep
        
      colnames(drop.prob) <- cell.names
      rownames(drop.prob) <- gene.names
      colnames(keep) <- cell.names
      rownames(keep) <- gene.names
        
      assays(sim)$DropProb <- drop.prob
      assays(sim)$Dropout <- !keep
        
      
      if(dir.exists(dirname)){
        if(verbose)cat("Saving counts with dropout shape ",dm, ": ", 100*round(mean(counts==0), 4),"% zeros...\n ")
        saveRDS(counts, paste0(dirname, trial, "dm", round(dm, 3), "_Raw.rds"))
      }
    }
    
    assays(sim)$counts = counts
    return(sim)
}

#' Get log-normal factors
#'
#' Randomly generate multiplication factors from a log-normal distribution.
#'
#' @param n.facs Number of factors to generate.
#' @param sel.prob Probability that a factor will be selected to be different
#'        from 1.
#' @param neg.prob Probability that a selected factor is less than one.
#' @param fac.loc Location parameter for the log-normal distribution.
#' @param fac.scale Scale factor for the log-normal distribution.
#'
#' @return Vector containing generated factors.
#' @importFrom stats rbinom rlnorm
getLNormFactors <- function(n.facs, is.selected, neg.prob=0, fac.loc, fac.scale) {

    n.selected <- sum(is.selected)
    if(length(fac.loc)>1){
      fac.loc = fac.loc[is.selected]
    }
    dir.selected <- (-1)^rbinom(n.selected, 1, neg.prob)
    if(length(fac.loc)>1){
      facs.selected <- sapply(1:n.selected, function(i)rlnorm(1, fac.loc[i], fac.scale))
    }
    else facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
    
    # Reverse directions for factors that are less than one
    dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected ^ dir.selected

    return(factors)
}


#' Simulate observed count matrix given technical biases and the true counts using methods in SymSim
#' @reference
#' Zhang X, Xu C, Yosef N. Simulating multiple faceted variability in single cell RNA sequencing. 
#' Nature communications. 2019 Jun 13;10(1):1-6. \url{https://www.nature.com/articles/s41467-019-10500-w}
#' 
#' @param true_counts gene cell matrix
#' @param meta_cell the meta information related to cells, will be combined with technical cell level information and returned 
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10 percent efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR1 the number of PCR cycles, default is 16
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @param hge2true if we add high gene expression to true counts
#' @param SE input, should be a summerized experiment rather than a list of elements, default is False
#' @param nbatch number of batches
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom utils data
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import progress
#' @rdname escoSimNoise
#' @export
Observed_Counts <- function(sim, trial, verbose, numCores =2, protocol = "UMI", nbatch=1){
  params <-metadata(sim)$Params
  dirname <- getParam(params, "dirname")
  true_counts = assays(sim)$TrueCounts
  nGenes <- getParam(params, "nGenes")
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, nGenes, replace = FALSE)
  rowData(sim)$genelen = gene_len
  
  alpha_mean_vec <- getParam(params, "alpha_mean")
  alpha_sd <- getParam(params, "alpha_sd")
  lenslope <- getParam(params, "lenslope")
  nbins <- getParam(params, "nbins")
  amp_bias_limit <- getParam(params, "amp_bias_limit")
  rate_2PCR <- getParam(params, "rate_2PCR")
  nPCR1 <- getParam(params, "nPCR1")
  nPCR2 <- getParam(params, "nPCR2") 
  LinearAmp <- getParam(params, "LinearAmp")
  LinearAmp_coef <- getParam(params, "LinearAmp_coef")
  depth_mean_vec <- getParam(params, "depth_mean")
  depth_sd <- getParam(params, "depth_sd")
  
  ngenes <- dim(true_counts)[1]
  ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005
  depth_lb <- 200 # lower bound for capture efficiency and sequencing depth  
  
  if((length(alpha_mean_vec)>1 | length(depth_mean_vec)>1) & (!dir.exists(dirname))){
    message("Warning: detect calls for simulating multiple configuration, but no directory to save files....")
    message("Modifying: only do simulation for the first configuration....")
    alpha_mean_vec = alpha_mean_vec[1]
    depth_mean_vec =  depth_mean_vec[1]
  }
  
  observed_counts = true_counts
  
  for(alpha_mean in alpha_mean_vec){
    for(depth_mean in depth_mean_vec){
      rate_2cap_vec <- rnorm_truc(n=ncells, mean = alpha_mean, sd=alpha_sd, a=rate_2cap_lb, b=Inf)
      depth_vec <- rnorm_truc(n=ncells, mean = depth_mean, sd=depth_sd,a=depth_lb, b=Inf)
      if(is.null(numCores))numCores=detectCores() -1
      cl <- makeCluster(numCores)
      registerDoSNOW(cl)
      total <- ncells
      pb <- progress_bar$new(
        format = "progress = :letter [:bar] :elapsed | eta: :eta", total = total, width = 60)
      progress <- function(n){
        pb$tick(tokens = list(letter = rep("", total)[n]))
      } 
      opts <- list(progress = progress)
      observed_counts <- foreach(i = c(1:total), .options.snow = opts, .export=c("amplify_cell")) %dopar% {
        return(amplify_cell(true_counts_1cell =  true_counts[, i], protocol=protocol, 
                     rate_2cap=rate_2cap_vec[i], gene_len=gene_len, amp_bias = amp_bias, 
                     rate_2PCR=rate_2PCR, nPCR1=nPCR1, nPCR2=nPCR2, LinearAmp = LinearAmp, 
                     LinearAmp_coef = LinearAmp_coef, N_molecules_SEQ = depth_vec[i]))     
      }
      stopCluster(cl)
      
      ## assign random batch ID to cells
      batchIDs <- sample(1:nbatch, ncells, replace = TRUE)
      meta_cell <- data.frame(alpha=rate_2cap_vec,depth=depth_vec, batch=batchIDs)
      
      if(protocol=="UMI"){
        UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
        nreads_perUMI <- lapply(observed_counts, "[[", 2)
        nUMI2seq <- sapply(observed_counts, "[[", 3)
        observed_counts <- UMI_counts
      } 
      else observed_counts <- do.call(cbind,observed_counts)
      
      # if (nbatch>1){
      #   mean_matrix <- matrix(0, ngenes, nbatch)
      #   batch_effect_size <- 2
      #   gene_mean <- rnorm(ngenes, 0, 1)
      #   temp <- lapply(1:ngenes, function(igene) {
      #     return(runif(nbatch, min = gene_mean[igene]-batch_effect_size, max = gene_mean[igene]+batch_effect_size))
      #   })
      #   mean_matrix <- do.call(rbind, temp)
      #   
      #   batch_factor <- matrix(0, ngenes, ncells)
      #   for (igene in 1:ngenes){
      #     for (icell in 1:ncells){
      #       batch_factor[igene, icell] <- rnorm(n=1, mean=mean_matrix[igene, batchIDs[icell]], sd=0.01)
      #     }
      #   }
      #   observed_counts <- 2^(log2(observed_counts)+batch_factor)
      # }
      
      rownames(observed_counts) = rownames(true_counts)
      colnames(observed_counts) = colnames(true_counts)
      
      if(dir.exists(dirname)){
        if(verbose)cat("Saving observed counts with alpha mean", round(alpha_mean, 3), " and depth mean", depth_mean, ":", 100*round(mean(observed_counts == 0), 4), "% zeros...\n")
        saveRDS(observed_counts, paste0(dirname, trial, "depme", depth_mean, "_alphame", round(alpha_mean,3), "_Raw.rds"))
      }
    }
  }
  assays(sim)$observedcounts = observed_counts
  return(sim)
  }


#' Simulate technical biases using methods in SymSim
#' @param lenslope amount of length bias. This value sould be less than 2*amp_bias_limit[2]/(nbins-1)
#' @param nbins number of bins for gene length
#' @param gene_len transcript length of each gene
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @export
cal_amp_bias <- function(lenslope, nbins, gene_len, amp_bias_limit){
  
  ngenes <- length(gene_len)
  len_bias_bin <- (-c(1:nbins))*lenslope
  len_bias_bin <- len_bias_bin-median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)
  
  rand_bias <- rnorm(ngenes, mean=0, sd=max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  #rand_bias <- runif(ngenes, -max_rand_bias,  max_rand_bias)
  
  binsize <- floor(ngenes/nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in 1:(nbins-1)){
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin-1)*binsize+1) : (ibin*binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins-1)*binsize+1) : ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins
  
  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias+len_bias
  return(amp_bias)
}


#' This function simulates the amplification, library prep, and the sequencing processes using methods in SymSim.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR1 the number of PCR cycles
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#' @importFrom utils data
#' @return read counts (if protocol="nonUMI") or UMI counts (if protocol="UMI)
#' @export
amplify_cell<- function(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias, 
                          rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ){
  
  # expand transcript counts to a vector of binaries of the same length of as the number of transcripts
  expandbinary<- function(true_counts_1cell){
    expanded_vec <- rep(1, sum(true_counts_1cell))
    trans_idx <- sapply(which(true_counts_1cell>0), 
                        function(igene){return(rep(igene, true_counts_1cell[igene]))})
    trans_idx <- unlist(trans_idx)
    return(list(expanded_vec, trans_idx))
  }
  
  ngenes <- length(gene_len)
  if (protocol=="nonUMI"){data(len2nfrag)} else 
    if(protocol=="UMI"){ } else
    {stop("protocol input should be nonUMI or UMI")}
  inds <- vector("list",2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expandbinary(c(true_counts_1cell,1))
  expanded_vec <- expanded_res[[1]]; trans_idx <- expanded_res[[2]]
  
  inds[[1]] <- which(expanded_vec > 0); expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]
  
  captured_vec <- expanded_vec; captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  if (sum(captured_vec[1:(length(captured_vec)-1)]) < 1) {return(rep(0, ngenes))}
  captured_vec[length(captured_vec)] <- 1
  inds[[2]] <- which(captured_vec > 0); captured_vec <- captured_vec[inds[[2]]]
  trans_idx <- trans_idx[inds[[2]]]
  
  amp_rate <- c((rate_2PCR+amp_bias[trans_idx[1:(length(trans_idx)-1)]]),1)
  
  # pre-amplification:
  if (LinearAmp){
    PCRed_vec <- captured_vec*LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    temp <- temp*2+captured_vec-temp
    for (iPCR in 2:nPCR1){
      eff <- runif(length(temp))*amp_rate
      v1 <- temp*(1-eff)
      round_down <- ((v1-floor(v1)) < runif(length(v1)))
      v1[round_down] <- floor(v1[round_down]); v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2*(temp-v1)
    }
    PCRed_vec <- temp
  }
  
  if (protocol=="nonUMI"){ # add fragmentation step here
    temp_vec <- PCRed_vec
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec; 
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    amp_mol_count=numeric(ngenes);
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]]
      amp_mol_count[i]=sum(x)
    }
    
    # for every copy of each transcript, convert it into number of fragments
    frag_vec <- numeric(ngenes)
    for (igene in which(amp_mol_count>0)){
      frag_vec[igene] <- sum(sample(len2nfrag[as.character(gene_len[igene]),], 
                                    amp_mol_count[igene], replace = TRUE))}
    # another 8 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    for (iPCR in 3:nPCR2){
      frag_vec <- frag_vec + round(frag_vec*rate_2PCR)
    }
    SEQ_efficiency=N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1) {read_count <- frag_vec} else{
      read_count <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)}) }
    return(read_count)
  } else if (protocol=="UMI"){
    
    
    get_prob <- function(glength){
      if (glength >= 1000){prob <- 0.7} else{
        if (glength >= 100 & glength < 1000){prob <- 0.78}
        else if (glength < 100) {prob <- 0}
      }
      return(prob)
    }
    prob_vec <- sapply(gene_len[trans_idx[1:(length(trans_idx)-1)]], get_prob)
    # fragmentation: 
    frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
    {return(rbinom(n=1, size = PCRed_vec[igene], prob = prob_vec[igene] ))})
    
    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    
    frag_vec <- round(frag_vec * (1+rate_2PCR)^(nPCR2-1))
    
    SEQ_efficiency <- N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1){sequenced_vec <- frag_vec} else {
      sequenced_vec <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)})}
    
    temp_vec <- c(sequenced_vec,1)
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec; 
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    
    UMI_counts=numeric(ngenes); 
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]];
      UMI_counts[i]=sum(x>0); 
    }
    
    return(list(UMI_counts, sequenced_vec, sum(frag_vec>0)))
  }
}

#' sample from truncated normal distribution (this function is borrowed from splatter)
#' @param a the minimum value allowed 
#' @param b the maximum value allowed
#' @export
rnorm_truc <- function(n, mean, sd, a, b){
  vec1 <- rnorm(n, mean = mean, sd=sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i){
      while (TRUE){
        temp <- rnorm(1, mean = mean, sd=sd)
        if (temp > a | temp > b) {break}}
      return(temp)} )
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}



#' Get Beta step probabilities
#'
#' Use a Beta distribution for set probabilities along a path (this function is borrowed from splatter)
#'
#' @param steps Number of steps
#' @param alpha Alpha parameter
#' @param beta Beta parameter
#' @importFrom S4Vectors metadata metadata<-
#'
#' @details
#' The density is sampled from a Beta distribution between 0 and 1. Infinite
#' densities at edges are adjusted and then the values are scaled to give
#' probabilities.
#'
#' @return Vector of probabilities
#'
#' @importFrom stats dbeta
#' @export
getBetaStepProbs <- function(steps, alpha, beta) {
  dens <- dbeta(seq(0, 1, length.out = steps), alpha, beta)
  
  # Adjust for infinite values at edge of distribution
  dens.inf <- !is.finite(dens)
  if (any(dens.inf) && all(dens[!dens.inf] == 0)) {
    dens[dens.inf] <- 1
  }
  if (!is.finite(dens[1])) {
    dens[1] <- 1.1 * dens[2]
  }
  if (!is.finite(dens[steps])) {
    dens[steps] <- 1.1 * dens[steps - 1]
  }
  
  probs <- dens / sum(dens)
  
  return(probs)
}

#' Sample density
#'
#' Sample from a density object using rejection sampling (this function is borrowed from splatter)
#'
#' @param n Number of values to sample
#' @param dens Density object to sample from
#' @param lower Lower x-axis bound on sampled values
#' @importFrom S4Vectors metadata metadata<-
#'
#' @details
#' Random points (x and y) are generated inside the range of the density object.
#' If they value is less than the density for that x value (and x is greater
#' than \code{lower}) then that x value is retained. Ten thousand points are
#' generated at a time until enough valid values have been sampled.
#'
#' @return Vector of sampled values
#'
#' @importFrom stats approxfun
#' @export
sampleDensity <- function(n, dens, lower = 0) {
  
  xmin <- min(dens$x)
  xmax <- max(dens$x)
  ymin <- min(dens$y)
  ymax <- max(dens$y)
  
  boundary <- approxfun(dens$x, dens$y)
  
  values <- c()
  nsel <- 0
  
  while(nsel < n) {
    x <- runif(1e4, xmin, xmax)
    y <- runif(1e4, ymin, ymax)
    sel <- y < boundary(x) & x > lower
    
    nsel <- nsel + sum(sel)
    values <- c(values, x[sel])
  }
  
  values <- values[seq_len(n)]
  
  return(values)
}


#' @export
randcop <-function(Rho, nCells){
  Col = chol(Rho)
  nGenes = nrow(Rho)
  copular = matrix(rnorm(nGenes*nCells), ncol = nCells)
  copular = t(Col) %*% copular
  copular = pnorm(copular)
  return(copular)
}

#' @export
makespd<-function(rho){
  er = eigen(rho)
  if(min(er$values)<0){
    oldsum = sum(er$values)
    er$values = er$values - min(er$values) + 1e-6
    newsum = sum(er$values)
    er$values = er$values/newsum*oldsum 
    rhocop = er$vectors %*% diag(er$values) %*% t(er$vectors)
    rhocop = rhocop%*%diag(1/diag(rhocop))
    return(rhocop)
  }
  else{
    return(rho)
  }
}


#' @importFrom utils data
#' @importFrom stats as.dist hclust cor
#' @export
randcor <- function(ngenes){
  data(puri_data)
  corr = cor(t(puri_data[sample(nrow(puri_data),ngenes),]))
  d <- stats::as.dist((1 - corr)/2)
  h <- stats::hclust(d)
  order <- h$order
  ans <- corr[order, order]
  ans = makespd(ans)
  return(ans)
}

#' Logistic function
#'
#' Implementation of the logistic function (this function is borrowed from splatter)
#'
#' @param x value to apply the function to.
#' @param x0 midpoint parameter. Gives the centre of the function.
#' @param k shape parameter. Gives the slope of the function.
#'
#' @return Value of logistic funciton with given parameters
logistic <- function(x, x0, k) {
  1 / (1 + exp(-k * (x - x0)))
}

#' Bind rows (matched)
#'
#' Bind the rows of two data frames, keeping only the columns that are common
#' to both (this function is borrowed from splatter).
#'
#' @param df1 first data.frame to bind.
#' @param df2 second data.frame to bind.
#'
#' @return data.frame containing rows from \code{df1} and \code{df2} but only
#'         common columns.
rbindMatched <- function(df1, df2) {
  common.names <- intersect(colnames(df1), colnames(df2))
  if (length(common.names) < 2) {
    stop("There must be at least two columns in common")
  }
  combined <- rbind(df1[, common.names], df2[, common.names])
  
  return(combined)
}

#' Winsorize vector
#'
#' Set outliers in a numeric vector to a specified percentile (this function is borrowed from splatter).
#'
#' @param x Numeric vector to winsorize
#' @param q Percentile to set from each end
#'
#' @return Winsorized numeric vector
winsorize <- function(x, q) {
  
  checkmate::check_numeric(x, any.missing = FALSE)
  checkmate::check_number(q, lower = 0, upper = 1)
  
  lohi <- stats::quantile(x, c(q, 1 - q), na.rm = TRUE)
  
  if (diff(lohi) < 0) { lohi <- rev(lohi) }
  
  x[!is.na(x) & x < lohi[1]] <- lohi[1]
  x[!is.na(x) & x > lohi[2]] <- lohi[2]
  
  return(x)
}

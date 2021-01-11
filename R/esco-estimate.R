#' Estimate esco simulation parameters
#'
#' Estimate simulation parameters for the esco simulation from a real
#' dataset. See the individual estimation functions for more details on how this
#' is done.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param dirname a tring of directory name to indicate where to save the results
#' @param group whether the data is believed to be of discrete cell groups or not, 
#'        if yes, the corresponding cellinfo indicating the cell group labels need to be input as well
#' @param cellinfo a vector of length n, where n is the number of cells. 
#'        Each entries is the group identity of a cell.
#' @param params escoParams object to store estimated values in.
#'
#' @seealso
#' \code{\link{escoEstMean}},  \code{\link{escoEstLib}},
#' \code{\link{escoEstOutlier}}, \code{\link{escoEstBCV}},
#' \code{\link{escoEstDropout}}, \code{\link{escoEstDE}},
#' \code{\link{escoEstGroupMean}}, \code{\link{escoEstGroupLib}},
#' \code{\link{escoEstGroupOutlier}}
#'
#' @return escoParams object containing the estimated parameters.
#'
#' @examples
#' # Load example data
#' library(scater)
#' data("sc_example_counts")
#'
#' params <- escoEstimate(sc_example_counts)
#' params
#' @rdname escoEstimate
#' @importFrom splatter setParams setParam getParams getParam 
#' @export
escoEstimate <- function(counts, dirname, group = FALSE, cellinfo = NULL, params = newescoParams()) {
    UseMethod("escoEstimate")
}

#' @rdname escoEstimate
#' @importFrom splatter setParams setParam getParams getParam
#' @export
#' 
escoEstimate.SingleCellExperiment <- function(counts, dirname, group = FALSE, cellinfo = NULL,
                                               params = newescoParams()) {
    counts <- BiocGenerics::counts(counts)
    escoEstimate(counts, params)
}


#' @rdname escoEstimate
#' @importFrom stats median quantile
#' @importFrom SC3 get_marker_genes
#' @importFrom splatter setParams setParam getParams getParam
#' @export
escoEstimate.matrix <- function(counts, dirname, group = FALSE, cellinfo = NULL, params = newescoParams()) {

  checkmate::assertClass(params, "escoParams")

  if(!group){
    # Normalise for library size and remove all zero genes
    lib.sizes <- colSums(counts)
    lib.medhigh <- median(lib.sizes[which(lib.sizes>quantile(lib.sizes, 0.1))])
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]
    norm.countshigh <- t(t(counts) / lib.sizes * lib.medhigh)
    norm.countshigh <- norm.countshigh[rowSums(norm.countshigh > 0) > 1, ]
    
    params <- escoEstDropout(norm.counts, params)
    
    params <- escoEstLib(counts, norm.counts, params)
    params <- escoEstMean(norm.counts, counts, params)
    params <- escoEstOutlier(norm.counts, params)
    params <- escoEstBCV(counts, norm.counts, params)
   
    params <- setParams(params, nGenes = nrow(counts),
                        nCells = ncol(counts))
    params <- setParams(params, dropout.type = "zeroinflate")
    return(params)
  }
  # if the cells are of discrete cell groups, the estimated differently
  else{
    lib.sizes <- colSums(counts)
    lib.med <- mean(lib.sizes)
    
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    saveRDS(norm.counts, paste0(dirname,  "normcounts.rds"))
    libadjust = colSums(counts)/colSums(norm.counts)
    
    ngroups = length(cellinfo)
    group.prob = as.vector(table(cellinfo))/ncol(counts)
    group.prob[1] = 1 - sum(group.prob[-1])
    params <- setParams(params, group.prob = group.prob)
   
    params <- escoEstGroupLib(counts = counts, params = params)
  
    markers = get_marker_genes(norm.counts, cellinfo)
    markers$genes = rownames(counts)
    save(markers, file = paste0(dirname, "markers.rdata"))
    degenes = markers[which(markers$auroc > 0.7), ]
    deall.prob <- nrow(degenes)/nrow(counts)
    de.prob <- as.vector(table(degenes$clusts))/nrow(degenes)
    params <- setParams(params, de.prob = de.prob, de.downProb = 0)
    
    housegenes <- markers$genes[which(markers$auroc < 0.2)]
    housegenes <- housegenes[which(rowMeans(norm.counts[housegenes,])> quantile(rowMeans(norm.counts), 0.8))]
    if(length(housegenes) > 0){
      house.prob = length(housegenes)/nrow(counts)
    }
    else{
      house.prob = 0.05
    }
  
    degenes.name = as.character(degenes$genes)
    nondegenes.name = setdiff(rownames(counts), as.character(degenes$genes))
    nonDE.normcounts = norm.counts[nondegenes.name,]
    nonDE.counts = counts[nondegenes.name,]
    DE.normcounts = norm.counts[degenes.name,]
    DE.counts = counts[degenes.name,]
  
    params <- escoEstGroupMean(nonDE.normcounts, DE.normcounts, degenes, cellinfo, params)
    params <- escoEstBCV(counts, norm.counts, params, cellinfo)
    params <- escoEstGroupOutlier(nonDE.normcounts, DE.normcounts, degenes, cellinfo, params)
    de <- escoEstDE(nonDE.normcounts, DE.normcounts, degenes, cellinfo, params)
    de.facLoc = de[[1]]
    de.facScale = de[[2]]
    de.rank = de[[3]]
    de.facrank = de[[4]]
    params <- setParams(params, de.facLoc = de.facLoc,
                        de.facScale  = de.facScale,
                        de.rank = de.rank,
                        de.facrank = de.facrank, nGroups = length(de.rank))
    params <- escoEstDropout(counts, params)
    params <- setParams(params, dropout.type = "zeroinflate")
    
    params <- setParams(params, nGenes = nrow(counts), nCells = ncol(counts), deall.prob = deall.prob, house.prob = house.prob)
    return(params)
  }
}

#' Estimate esco library size parameters
#'
#' The Shapiro-Wilks test is used to determine if the library sizes are
#' normally distributed. If so a normal distribution is fitted to the library
#' sizes, if not (most cases) a log-normal distribution is fitted and the
#' estimated parameters are added to the params object. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param norm.counts library size normalised counts matrix.
#' @param params escoParams object to store estimated values in.
#'
#' @return escoParams object with estimated values.
#' @rdname escoEstLib
#' @importFrom stats shapiro.test
#' @importFrom splatter setParams setParam getParams getParam
escoEstLib <- function(counts, norm.counts, params) {
  
  #lib.sizes <- colSums(counts)
  drop.mid <- getParam(params, "dropout.mid")
  drop.shape <- getParam(params, "dropout.shape")
  
  #means <- winsorize(means, q = 0.01)
  means = rowMeans(norm.counts)
  means = means[which(means!=0)]
  lmeans <- log(means)
  drop.prob = logistic(lmeans, x0 = drop.mid, k = drop.shape)
  #keep.prob = (1 - drop.prob)
  keep.prob = (1 - drop.prob)/(1-dpois(0, lambda = means))
  keep.prob[keep.prob<0] = 1
  keep.prob[is.na(keep.prob)] = 1
  keep.prob[keep.prob>1] = 1
  
  truecounts = counts/keep.prob
  
  lib.sizes = colSums(truecounts)
    
  #lib.sizes <- lib.sizes[which(lib.sizes>quantile(lib.sizes, 0.1))]
  
  
  if (length(lib.sizes) > 5000) {
    message("NOTE: More than 5000 cells provided. ",
            "5000 sampled library sizes will be used to test normality.")
    lib.sizes.sampled <- sample(lib.sizes, 5000, replace = FALSE)
  } else {
    lib.sizes.sampled <- lib.sizes
  }
  
  norm.test <- shapiro.test(lib.sizes.sampled)
  lib.norm <- norm.test$p.value > 0.2
  
  
  if (lib.norm) {
    fit <- fitdistrplus::fitdist(lib.sizes, "norm")
    lib.loc <- unname(fit$estimate["mean"])
    lib.scale <- unname(fit$estimate["sd"])
    message("NOTE: Library sizes have been found to be normally ",
            "distributed instead of log-normal. You may want to check ",
            "this is correct.")
  } else {
    fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")
    lib.loc <- unname(fit$estimate["meanlog"])
    lib.scale <- unname(fit$estimate["sdlog"])
  }
  # 
  # if(gamma){
  #   fit <- fitdistrplus::fitdist(lib.sizes, "gamma", method = "mge",
  #                                gof = "CvM")
  #   if (fit$convergence > 0) {
  #     warning("Fitting means using the Goodness of Fit method failed, ",
  #             "using the Method of Moments instead")
  #     fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")
  #   }
  #   lib.loc = unname(fit$estimate["shape"])
  #   lib.scale = unname(fit$estimate["rate"])
  # }
  
  params <- setParams(params, lib.loc = lib.loc, lib.scale = lib.scale,
                      lib.norm = lib.norm, lib.dens = density(lib.sizes))
  
  return(params)
}


#' Estimate esco group mean parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate gene expression means.
#'
#' @param normcounts library size normalised counts matrix.
#' @param counts counts matrix to estimate parameters from.
#' @param params escoParams object to store estimated values in.
#'
#' @details
#' Parameter for the gamma distribution are estimated by fitting the mean
#' normalised counts using \code{\link[fitdistrplus]{fitdist}}. The 'maximum
#' goodness-of-fit estimation' method is used to minimise the Cramer-von Mises
#' distance. This can fail in some situations, in which case the 'method of
#' moments estimation' method is used instead. Prior to fitting the means are
#' winsorized by setting the top and bottom 10 percent of values to the 10th
#' and 90th percentiles.
#' @import DescTools
#' @importFrom splatter setParams setParam getParams getParam
#' @rdname escoEstMean
#' @return escoParams object with estimated values.
escoEstMean <- function(normcounts, counts, params) {
  
  means <- rowMeans(normcounts)
  means <- means[means != 0]
  
  drop.mid <- getParam(params, "dropout.mid")[1]
  drop.shape <- getParam(params, "dropout.shape")
  
  lmeans <- log(means)
  #means <- winsorize(means, q = 0.01)
  
  drop.prob = logistic(lmeans, x0 = drop.mid, k = drop.shape)
  #keep.prob = (1 - drop.prob)
  keep.prob = (1 - drop.prob)/(1-dpois(0, lambda = means))
  keep.prob[keep.prob<0] = 0
  keep.prob[is.na(keep.prob)] = 0
  keep.prob[keep.prob>1] = 1
  means = means/keep.prob
  
  lmeans <- log(means)
  
  #means <- means[which(means>quantile(means, prob = 0.1))]
  
  fit <- fitdistrplus::fitdist(means, "gamma", method = "mge",
                               gof = "CvM")
  if (fit$convergence > 0) {
    warning("Fitting means using the Goodness of Fit method failed, ",
            "using the Method of Moments instead")
    fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")
  }
  
  params <- setParams(params, mean.shape = unname(fit$estimate["shape"]),
                      mean.rate = unname(fit$estimate["rate"]), mean.dens = density(lmeans))
  
  return(params)
}



#' Estimate esco expression outlier parameters for a single cell group
#'
#' Parameters are estimated by comparing means of individual genes to the
#' median mean expression level.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params escoParams object to store estimated values in.
#'
#' @details
#' Expression outlier genes are detected using the Median Absolute Deviation
#' (MAD) from median method. If the log2 mean expression of a gene is greater
#' than two MADs above the median log2 mean expression it is designated as an
#' outlier. The proportion of outlier genes is used to estimate the outlier
#' probability. Factors for each outlier gene are calculated by dividing mean
#' expression by the median mean expression. A log-normal distribution is then
#' fitted to these factors in order to estimate the outlier factor location and
#' scale parameters using \code{\link[fitdistrplus]{fitdist}}.
#' @rdname escoEstOutlier
#' @return escoParams object with estimated values.
#' @importFrom splatter setParams setParam getParams getParam
escoEstOutlier <- function(norm.counts, params) {
  
  means <- rowMeans(norm.counts)
  lmeans <- log(means)
  
  med <- median(lmeans)
  mad <- mad(lmeans)
  
  bound <- med + 2*mad
  
  outs <- which(lmeans > bound)
  
  prob <- length(outs) / nrow(norm.counts)
  
  params <- setParams(params, out.prob = prob)
  
  if (length(outs) > 1) {
    facs <- means[outs] / median(means)
    fit <- fitdistrplus::fitdist(facs, "lnorm")
    
    params <- setParams(params,
                        out.facLoc = unname(fit$estimate["meanlog"]),
                        out.facScale = unname(fit$estimate["sdlog"]))
  }
  
  return(params)
}

#' Estimate esco Biological Coefficient of Variation parameters
#'
#' Parameters are estimated using the \code{\link[edgeR]{estimateDisp}} function
#' in the \code{edgeR} package.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param norm.counts normalized counts matrix to estimate parameters from.
#' @param params escoParams object to store estimated values in.
#' @param cellinfo info about the identity of each cell in the cell structure.
#'        If cellinfo is not null, then BCV is estiamted adjusted to the cell structures 
#' @details
#' The \code{\link[edgeR]{estimateDisp}} function is used to estimate the common
#' dispersion and prior degrees of freedom. See
#' \code{\link[edgeR]{estimateDisp}} for details. When estimating parameters on
#' simulated data we found a broadly linear relationship between the true
#' underlying common dispersion and the \code{edgR} estimate, therefore we
#' apply a small correction, \code{disp = 0.1 + 0.25 * edgeR.disp}.
#' @rdname escoEstBCV
#' @return escoParams object with estimated values.
#' @importFrom splatter setParams setParam getParams getParam
escoEstBCV <- function(counts, norm.counts, params, cellinfo = NULL){
  if(is.null(cellinfo)){
    means = rowMeans(norm.counts)
    drop.mid <- getParam(params, "dropout.mid")[1]
    drop.shape <- getParam(params, "dropout.shape")
    
    #means <- winsorize(means, q = 0.01)
    lmeans <- log(means)
    drop.prob = logistic(lmeans, x0 = drop.mid, k = drop.shape)
    keep.prob = (1 - drop.prob)
    #keep.prob = (1 - drop.prob)/(1-dpois(0, lambda = means))
    keep.prob[keep.prob<0] = 0
    keep.prob[is.na(keep.prob)] = 0
    keep.prob[keep.prob>1] = 1
    counts = counts/keep.prob
  }
  
  if(is.null(cellinfo)){
    design <- matrix(1, ncol(counts), 1) 
    disps <- edgeR::estimateDisp(counts, design = design)
    params <- setParams(params,
                        bcv.common = 0.1 + 0.25 * disps$common.dispersion,
                        bcv.df = disps$prior.df)
  }
  else{
    disps <- edgeR::estimateDisp(counts, group = cellinfo)
    params <- setParams(params,
                        bcv.common = 0.1 + 0.25 * disps$common.dispersion,
                        bcv.df = disps$prior.df)
  }
  
  
  return(params)
}


#' Estimate esco dropout parameters
#'
#' Estimate the midpoint and shape parameters for the logistic function used
#' when simulating dropout.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params escoParams object to store estimated values in.
#'
#' @details
#' Logistic function parameters are estimated by fitting a logistic function
#' to the relationship between log2 mean gene expression and the proportion of
#' zeros in each gene. See \code{\link[stats]{nls}} for details of fitting.
#' Note this is done on the experiment level.
#'
#' @return escoParams object with estimated values.
#' @rdname escoEstDropout
#' @importFrom stats dnbinom nls
#' @importFrom splatter setParams setParam getParams getParam
escoEstDropout <- function(norm.counts, params) {
  
  means <- rowMeans(norm.counts)
  
  x <- log(means)
  
  obs.zeros <- rowSums(norm.counts == 0)
  
  y <- obs.zeros / ncol(norm.counts)
  
  df <- data.frame(x, y)
  
  x_approx_mid <- median(x[which(y > 0.0001 & y < 0.9999)])
  fit <- nls(y ~ logistic(x, x0 = x0, k = k), data = df,
             start = list(x0 = x_approx_mid, k = -1))
  
  mid <- summary(fit)$coefficients["x0", "Estimate"]
  shape <- summary(fit)$coefficients["k", "Estimate"]
  
  params <- setParams(params, dropout.mid = mid, dropout.shape = shape)
  
  return(params)
}




#' Estimate esco mean parameters for discrete cell groups
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate gene expression means for each cell group.
#'
#' @param nonDE.normcounts library size normalised counts matrix for nonDE genes.
#' @param DE.normcounts library size normalised counts matrix for DE genes.
#' @param degenes a dataframe of two columns: "genes" and "clust", 
#'        where "genes" contains the position of a degenes, and "clust" contains the 
#'        corresponding cell group the degenes marks.
#' @param cellinfo a vector of length n, where n is the number of cells. 
#'        Each entries is the group identity of a cell.
#' @param params escoParams object to store estimated values in.
#' 
#' @details
#' Parameter for the gamma distribution are estimated by fitting the mean
#' normalised counts using \code{\link[fitdistrplus]{fitdist}}. Particularly, 
#' the entries corresponds to DE genes in the cell group it marks are not considered in the fitting. 
#' The 'maximum goodness-of-fit estimation' method is used to minimise the Cramer-von Mises
#' distance. This can fail in some situations, in which case the 'method of
#' moments estimation' method is used instead. Prior to fitting the means are
#' winsorized by setting the top and bottom 10 percent of values to the 10th
#' and 90th percentiles.
#' @rdname escoEstGroupMean
#' @return escoParams object with estimated values.
#' @importFrom splatter setParams setParam getParams getParam
escoEstGroupMean <- function(nonDE.normcounts, DE.normcounts, degenes, cellinfo, params) {
    
    means <- rowMeans(nonDE.normcounts)
    for(idx in seq_len(unique(degenes$clusts))){
      marker = degenes$genes[which(degenes$clusts==idx)]
      cells = colnames(DE.normcounts)[which(cellinfo!=levels(cellinfo)[idx])]
      means = c(means, rowMeans(DE.normcounts[marker,cells]))
    }
    
    means <- means[means != 0]

    means <- winsorize(means, q = 0.01)
    lmeans <- log(means)

    fit <- fitdistrplus::fitdist(means, "gamma", method = "mge",
                                 gof = "CvM")
    if (fit$convergence > 0) {
        warning("Fitting means using the Goodness of Fit method failed, ",
                "using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")
    }

    params <- setParams(params, mean.shape = unname(fit$estimate["shape"]),
                        mean.rate = unname(fit$estimate["rate"]), mean.dens = density(lmeans))

    return(params)
}


#' Estimate esco library size parameters for discrete cell groups
#'
#' The Shapiro-Wilks test is used to determine if the library sizes are
#' normally distributed. If so a normal distribution is fitted to the library
#' sizes, if not (most cases) a log-normal distribution is fitted and the
#' estimated parameters are added to the params object. Specifically, for a cell group, the fitting
#' process uses only the cell samples within this cell group. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param cellinfo a vector of length n, where n is the number of cells. 
#'        Each entries is the group identity of a cell.
#' @param params escoParams object to store estimated values in.
#' 
#' @return escoParams object with estimated values.
#' @rdname escoEstGroupLib
#' @importFrom stats shapiro.test
#' @importFrom splatter setParams setParam getParams getParam
escoEstGroupLib <- function(counts, cellinfo = NULL, params) {
  
  if(!is.null(cellinfo)){
    lib.norm = c()
    lib.loc = c()
    lib.scale = c()
    for(idx in seq_len(length(unique(cellinfo)))){
      cells = colnames(counts)[which(cellinfo==levels(cellinfo)[idx])]
      lib.sizes <- colSums(counts[,cells])
      if (length(lib.sizes) > 5000) {
        message("NOTE: More than 5000 cells provided. ",
                "5000 sampled library sizes will be used to test normality.")
        lib.sizes.sampled <- sample(lib.sizes, 5000, replace = FALSE)
      } else {
        lib.sizes.sampled <- lib.sizes
      }
      
      norm.test <- shapiro.test(lib.sizes.sampled)
      lib.norm[idx] <- norm.test$p.value > 0.2
      
      if (lib.norm[idx]) {
        fit <- fitdistrplus::fitdist(lib.sizes, "norm")
        lib.loc[idx] <- unname(fit$estimate["mean"])
        lib.scale[idx] <- unname(fit$estimate["sd"])
        message("NOTE: Library sizes have been found to be normally ",
                "distributed instead of log-normal. You may want to check ",
                "this is correct.")
      } else {
        fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")
        lib.loc[idx] <- unname(fit$estimate["meanlog"])
        lib.scale[idx] <- unname(fit$estimate["sdlog"])
      }
    }
    
  }
  else{
    lib.sizes <- colSums(counts)
    if (length(lib.sizes) > 5000) {
      message("NOTE: More than 5000 cells provided. ",
              "5000 sampled library sizes will be used to test normality.")
      lib.sizes.sampled <- sample(lib.sizes, 5000, replace = FALSE)
    } else {
      lib.sizes.sampled <- lib.sizes
    }
    
    norm.test <- shapiro.test(lib.sizes.sampled)
    lib.norm <- norm.test$p.value > 0.2
    
    if (lib.norm) {
      fit <- fitdistrplus::fitdist(lib.sizes, "norm")
      lib.loc <- unname(fit$estimate["mean"])
      lib.scale <- unname(fit$estimate["sd"])
      message("NOTE: Library sizes have been found to be normally ",
              "distributed instead of log-normal. You may want to check ",
              "this is correct.")
    } else {
      fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")
      lib.loc <- unname(fit$estimate["meanlog"])
      lib.scale <- unname(fit$estimate["sdlog"])
    }
  }
  
  params <- setParams(params, lib.loc = lib.loc, lib.scale = lib.scale,
                      lib.norm = lib.norm)
  
  return(params)
}



#' Estimate esco expression DE parameters for discrete cell groups
#'
#' Parameters are estimated by comparing means of individual genes to the
#' median mean expression level.
#'
#' @param nonDE.normcounts library size normalised counts matrix for nonDE genes.
#' @param DE.normcounts library size normalised counts matrix for DE genes.
#' @param degenes a dataframe of two columns: "genes" and "clust", 
#'        where "genes" contains the position of a degenes, and "clust" contains the 
#'        corresponding cell group the degenes marks.
#' @param cellinfo a vector of length n, where n is the number of cells. 
#'        Each entries is the group identity of a cell.
#' @param params escoParams object to store estimated values in.
#' @rdname escoEstDE
#' @details
#' Expression DE factor are detected using the Median Absolute Deviation
#' (MAD) from median method. If the log2 mean expression of a gene is greater
#' than two MADs above the median log2 mean expression it is designated as an
#' outlier. The proportion of outlier genes is used to estimate the outlier
#' probability. Factors for each outlier gene are calculated by dividing mean
#' expression by the median mean expression. A log-normal distribution is then
#' fitted to these factors in order to estimate the outlier factor location and
#' scale parameters using \code{\link[fitdistrplus]{fitdist}}.
#'
#' @return escoParams object with estimated values.
#' @importFrom splatter setParams setParam getParams getParam
escoEstDE <- function(nonDE.normcounts, DE.normcounts, degenes, cellinfo, params) {
  de.rank = list()
  de.facrank = list()
  means <- rowMeans(nonDE.normcounts)
  for(idx in seq_len(length(unique(degenes$clusts)))){
    marker = degenes$genes[which(degenes$clusts==idx)]
    cells = colnames(DE.normcounts)[which(cellinfo!=levels(cellinfo)[idx])]
    means = c(means, rowMeans(DE.normcounts[marker,cells]))
  }
  
  ordall = sort(means, index.return=TRUE)$ix
  
  de.facLoc = c()
  de.facScale = c()
  
  start = nrow(nonDE.normcounts)+1
  
  for(idx in seq_len(length(unique(degenes$clusts)))){
    marker = as.character(degenes$genes[which(degenes$clusts==idx)])
    nonDEcells = colnames(DE.normcounts)[which(cellinfo!=levels(cellinfo)[idx])]
    nonDEmeans = rowMeans(DE.normcounts[marker, nonDEcells])
    
    DEcells = colnames(DE.normcounts)[which(cellinfo==levels(cellinfo)[idx])]
    Demeans = rowMeans(DE.normcounts[marker, DEcells])
    
    de.rank[[idx]] = ordall[start: (start + length(marker)-1)]
    de.facrank[[idx]] = sort(Demeans / nonDEmeans, index.return = TRUE)$ix
    start = start + length(marker)
    
    med <- mean(nonDEmeans)
    
    facs <- Demeans / med
    fit <- fitdistrplus::fitdist(facs, "lnorm")
    de.facLoc[idx] = unname(fit$estimate["meanlog"])
    de.facScale[idx] = unname(fit$estimate["sdlog"])
  }
  
  return(list(de.facLoc = de.facLoc, de.facScale = de.facScale, de.rank = de.rank, de.facrank = de.facrank))
}


#' Estimate esco Group expression outlier parameters for discrete cell groups
#'
#' Parameters are estimated by comparing means of individual genes to the
#' median mean expression level.
#'
#' @param nonDE.normcounts library size normalised counts matrix for non DE genes.
#' @param DE.normcounts library size normalised counts matrix for DE genes.
#' @param degenes a dataframe of two columns: "genes" and "clust", 
#'        where "genes" contains the position of a degenes, and "clust" contains the 
#'        corresponding cell group the degenes marks.
#' @param cellinfo a vector of length n, where n is the number of cells. 
#'        Each entries is the group identity of a cell.
#' @param params escoParams object to store estimated values in.
#' @rdname escoEstGroupOutlier
#' @details
#' Expression outlier genes are detected using the Median Absolute Deviation
#' (MAD) from median method. If the log2 mean expression of a gene is greater
#' than two MADs above the median log2 mean expression it is designated as an
#' outlier. Specifically, the meadian log2 mean expression for a degene is computed 
#' with the cell samples that does not belongs to the cell group it marked. 
#' The proportion of outlier genes is used to estimate the outlier
#' probability. Factors for each outlier gene are calculated by dividing mean
#' expression by the median mean expression. A log-normal distribution is then
#' fitted to these factors in order to estimate the outlier factor location and
#' scale parameters using \code{\link[fitdistrplus]{fitdist}}.
#'
#' @return escoParams object with estimated values.
#' @importFrom splatter setParams setParam getParams getParam
escoEstGroupOutlier <- function(nonDE.normcounts, DE.normcounts, degenes, cellinfo, params) {
    means <- rowMeans(nonDE.normcounts)
    for(idx in seq_len(length(unique(degenes$clusts)))){
      marker = degenes$genes[which(degenes$clusts==idx)]
      cells = colnames(DE.normcounts)[which(cellinfo!=levels(cellinfo)[idx])]
      means = c(means, rowMeans(DE.normcounts[marker,cells]))
    }
    #means <- rowMeans(norm.counts)
    lmeans <- log(means)

    med <- median(lmeans)
    mad <- mad(lmeans)

    bound <- med + 2 * mad

    outs <- which(lmeans > bound)

    prob <- length(outs) / (nrow(nonDE.normcounts) + nrow(DE.normcounts))

    params <- setParams(params, out.prob = prob)

    if (length(outs) > 1) {
        facs <- means[outs] / median(means)
        fit <- fitdistrplus::fitdist(facs, "lnorm")

        params <- setParams(params,
                            out.facLoc = unname(fit$estimate["meanlog"]),
                            out.facScale = unname(fit$estimate["sdlog"]))
    }

    return(params)
}

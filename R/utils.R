#' Sample from truncated normal distribution 
#' 
#' This function is borrowed from \code{\link[splatter]{splatter}}.
#' 
#' @param a the minimum value allowed 
#' @param b the maximum value allowed
#' @param mean the mean of the normal distribution
#' @param sd the standard deviation of the normal distribution
#' @param n sample size
#' @return a vector of size n 
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
#' Use a Beta distribution for set probabilities along a path. 
#' (this function is borrowed from \code{\link[splatter]{splatter}}).
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
#' Sample from a density object using rejection sampling 
#' (this function is borrowed from \code{\link[splatter]{splatter}}).
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

#' Sample from a gaussian copula 
#' @param Rho the correlation matrix in the gaussian copula
#' @param nCells the number of samples to draw
#' @return matrix of nCell rows, where each row is a sample from the gaussian copula
randcop <-function(Rho, nCells){
  Col = chol(Rho)
  nGenes = nrow(Rho)
  copular = matrix(rnorm(nGenes*nCells), ncol = nCells)
  copular = t(Col) %*% copular
  copular = pnorm(copular)
  return(copular)
}

#' Make positive definite 
#' Make a matrix positive definite  through tweaking the eigen values
#' @param rho the input matrix (assuming squared matrix) 
#' @return a positive definite matrix 
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

#' Random sampeled correlation from the reallife gene correlations
#'
#' Generate a random correlation matrix that is positive definite, using a purified gene expression data set 
#' \url{https://www.eurekalert.org/pub_releases/2017-11/sfn-nwa111417.php}. 
#' 
#' @importFrom utils data
#' @importFrom stats as.dist hclust cor
#' @param ngenes  the number of genes 
#' @return a correlation matrix of dimension ngenes*ngenes
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
#' Implementation of the logistic function  
#' (this function is borrowed from \code{\link[splatter]{splatter}}).
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
#' to both (this function is borrowed from \code{\link[splatter]{splatter}}).
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
#' Set outliers in a numeric vector to a specified percentile 
#' (this function is borrowed from \code{\link[splatter]{splatter}}).
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

#' Get log-normal factors
#'
#' Randomly generate multiplication factors from a log-normal distribution 
#' (this function is borrowed from \code{\link[splatter]{splatter}}).
#'
#' @param n.facs Number of factors to generate.
#' @param is.selected whether a factor will be selected to be different
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

maxfunc <- function(x, data){
  if(is.null(x))return(max(data))
  else return(x)
}

minfunc <- function(x, data){
  if(is.null(x))return(min(data))
  else return(x)
}

# library("copula")
# library("stats")
# library("VineCopula")
# 
# estimate_copula<-function(data) { 
#   corr_mat <- matrix(1, nrow= nrow(data), ncol = ncol(data))
#   data_col <- matrix(c(0,0), nrow=ncol(data), ncol = 2)
#   for (i in 1:(nrow(data)-1)){
#     for (j in (i+1):nrow(data)){
#       data_col[,1]=t(data[i,])
#       data_col[,2]=t(data[j,])
#       m = pobs(data_col)
#       cop = normalCopula(dim = 2)
#       fit = fitCopula(cop, m, method = 'ml')
#       corr_mat[i,j] = coef(fit)[1]
#       corr_mat[j,i] = corr_mat[i,j] 
#     }
#   }
#   return(corr_mat)
# }

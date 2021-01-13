# This file contains all the objected orientated related definitions

#====== Classes
#' The Params virtual class
#'
#' Virtual S4 class that all other Params classes inherit from.
#' @section Parameters:
#'
#' The Params class defines the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data.
#' @name Params
#' @rdname Params
#' @aliases Params-class
setClass("Params",
        contains = "VIRTUAL",
        slots = c(nGenes = "numeric",
                   nCells = "numeric",
                   seed = "numeric"),
        prototype = prototype(nGenes = 10000, nCells = 100,
                              seed = sample(1:1e6, 1)))

#' The escoParams class
#'
#' S4 class that holds parameters for the ESCO simulation.
#'
#' @section Parameters: 
#' The ESCO simulation requires the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Mean parameters}}{
#'         \describe{
#'             \item{\code{[mean.method]}}{Whether to use a gamma distribution 
#'             or a given density.}
#'             \item{\code{mean.dens}}{A density object.}
#'             \item{\code{mean.shape}}{Shape parameter for the mean gamma
#'             distribution.}
#'             \item{\code{mean.rate}}{Rate parameter for the mean gamma
#'             distribution.}
#'         }
#'     }
#'     \item{\emph{Library size parameters}}{
#'         \describe{
#'             \item{[lib.method]}{Whether to use a gamma distribution 
#'             or a given density.}
#'             \item{lib.dens}{A density object.}
#'             \item{\code{lib.loc}}{Location (meanlog) parameter for the
#'             library size log-normal distribution, or mean parameter if a
#'             normal distribution is used.}
#'             \item{\code{lib.scale}}{Scale (sdlog) parameter for the library
#'             size log-normal distribution, or sd parameter if a normal
#'             distribution is used.}
#'             \item{\code{lib.norm}}{Logical. Whether to use a normal
#'             distribution for library sizes instead of a log-normal.}
#'         }
#'     }
#'     \item{\emph{Expression outlier parameters}}{
#'         \describe{
#'             \item{\code{out.prob}}{Probability that a gene is an expression
#'             outlier.}
#'             \item{\code{out.facLoc}}{Location (meanlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'             \item{\code{out.facScale}}{Scale (sdlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'         }
#'     }
#'     \item{\emph{Group parameters}}{
#'         \describe{
#'             \item{\code{[nGroups]}}{The number of groups to simulate.}
#'             \item{\code{[group.prob]}}{Probability that a cell comes from a
#'             group.}
#'         }
#'     }
#'    
#'    \item{\emph{Tree parameters}}{
#'         \describe{
#'             \item{\code{[tree]}}{The tree structure to simulate.}
#'         }
#'     }
#'     
#'     \item{\emph{Differential expression parameters}}{
#'         \describe{
#'             \item{\code{[de.center]}}{The mean of the tree DE factors.}
#'             \item{\code{[de.prob]}}{Probability that a gene is differentially
#'             expressed in a group. Can be a vector.}
#'             \item{\code{[de.loProb]}}{Probability that a differentially
#'             expressed gene is down-regulated. Can be a vector.}
#'             \item{\code{[de.facLoc]}}{Location (meanlog) parameter for the
#'             differential expression factor log-normal distribution. Can be a
#'             vector.}
#'             \item{\code{[de.facScale]}}{Scale (sdlog) parameter for the
#'             differential expression factor log-normal distribution. Can be a
#'             vector.}
#'         }
#'     }
#'     \item{\emph{Biological Coefficient of Variation parameters}}{
#'         \describe{
#'             \item{\code{bcv.common}}{Underlying common dispersion across all
#'             genes.}
#'             \item{\code{bcv.df}}{Degrees of Freedom for the BCV inverse
#'             chi-squared distribution.}
#'         }
#'     }
#'     \item{\emph{Dropout parameters}}{
#'         \describe{
#'             \item{\code{[dropout.type]}}{The type of dropout to simulate.
#'             "none" indicates no dropout, "zeroinflate" uses zero 
#'             inflation model to add dropouts,
#'              "downsampling" uses similar procedure in 
#'              SymSim to mimic the experimental steps 
#'              for adding dropouts.}
#'             \item{\code{dropout.mid}}{Midpoint parameter for the dropout
#'             logistic function.}
#'             \item{\code{dropout.shape}}{Shape parameter for the dropout
#'             logistic function.}
#'             \item{\code{[alpha_mean]}}{Mean parameter 
#'             for the dwonsampling
#'              gamma function.}
#'             \item{\code{[alpha_sd]}}{Standard variance parameter 
#'             for the downsampling gamma function.}
#'             \item{\code{[lenslope]}}{Shape parameter for the dropout
#'             logistic function.}
#'              \item{\code{[nbins]}}{Shape parameter for the dropout
#'             logistic function.}
#'             \item{\code{[amp_bias_limt]}}{Shape parameter for the dropout
#'             logistic function.}
#'             \item{\code{[rate_2PCR]}}{PCR efficiency, usually very high}
#'             \item{\code{[LinearAmp]}}{if linear amplification is 
#'             used for pre-amplification step, default is FALSE}
#'             \item{\code{[LinearAmp_coef]}}{the coeficient of
#'              linear amplification, 
#'             that is, how many times each molecule is amplified by}
#'             \item{\code{[depth_mean]}}{Mean parameter of the 
#'             sequencing depths.}
#'             \item{\code{[depth_sd]}}{Standard variance parameter of 
#'             sequencing depths.}
#'         }
#'     }
#'  }
#'
#' The parameters not shown in brackets can be 
#' estimated from real data using
#' \code{\link{escoEstimate}}. For details of 
#' the Splatter simulation
#' see \code{\link{escoSimulate}}.
#'
#' @name escoParams
#' @rdname escoParams
#' @aliases escoParams-class
#' @exportClass escoParams
setClass("escoParams",
         contains = "Params",
         slots = c(mean.shape = "numeric",
                   mean.rate = "numeric",
                   mean.dens = "density",
                   mean.method = "character",
                   lib.loc = "numeric",
                   lib.scale = "numeric",
                   lib.norm = "logical",
                   lib.dens = "density",
                   lib.method = "character",
                   out.prob = "numeric",
                   out.facLoc = "numeric",
                   out.facScale = "numeric",
                   nGroups = "numeric",
                   group.prob = "numeric",
                   deall.prob = "numeric",
                   de.prob = "numeric",
                   de.quantile = "numeric",
                   de.downProb = "numeric",
                   de.facLoc = "numeric",
                   de.facScale = "numeric",
                   de.rank = "list",
                   de.facrank = "list",
                   paths.deprob = "numeric",
                   paths.design = "data.frame",
                   paths.factors = "list",
                   paths.means = "list",
                   paths.DEgenes = "numeric",
                   cells.design = "data.frame",
                   house.prob = "numeric",
                   bcv.common = "numeric",
                   bcv.df = "numeric",
                   dropout.type = "character",
                   dropout.mid = "numeric",
                   dropout.shape = "numeric",
                   dropout.cort = "logical",
                   withcorr = "logical",
                   corr = "list",
                   corr.prob = "numeric",
                   tree = "list",
                   de.center = "numeric",
                   alpha_mean = "numeric", 
                   alpha_sd = "numeric",
                   lenslope = "numeric", 
                   nbins= "numeric", 
                   amp_bias_limit =  "numeric",
                   rate_2PCR = "numeric", 
                   nPCR1 =  "numeric", 
                   nPCR2 =  "numeric", 
                   LinearAmp = "logical", 
                   LinearAmp_coef = "numeric", 
                   depth_mean =  "numeric", 
                   depth_sd =  "numeric",
                   dirname = "character",
                   trials = "numeric"
         ),
         prototype = prototype(mean.rate = 0.3,
                               mean.shape = 0.6,
                               mean.dens = density(log(rgamma(10, rate = 0.3,
                                                              shape = 0.6))),
                               mean.method = "fit",
                               lib.loc = 9,
                               lib.scale = 0.2,
                               lib.norm = FALSE,
                               lib.dens = density(rlnorm(10, 9, 0.2)),
                               lib.method = "fit",
                               out.prob = 0.05,
                               out.facLoc = 4,
                               out.facScale = 0.5,
                               nGroups = 1,
                               group.prob = 1,
                               deall.prob = 0.3,
                               de.prob = 0.1,
                               de.quantile = 0.7,
                               de.downProb = 0.5,
                               de.facLoc = 0.1,
                               de.facScale = 0.7,
                               de.rank = list(),
                               de.facrank = list(),
                               paths.deprob = 0.2,
                               paths.design = data.frame(
                                 Path = c(1, 2, 3),
                                 From = c(0, 1, 1),
                                 Steps = c(100, 100, 100)
                               ),
                               paths.factors = list(),
                               paths.DEgenes = c(0),
                               paths.means = list(),
                               cells.design = data.frame(
                                 Path = c(1,2,3),
                                 Probability = c(0.3, 0.3, 0.4),
                                 Alpha = 1,
                                 Beta = 1
                               ),
                               house.prob = 0.1,
                               bcv.common = 0.1,
                               bcv.df = 60,
                               dropout.type = c("downsample", "zeroinflate"),
                               dropout.mid = 1,
                               dropout.shape = -1,
                               dropout.cort = FALSE,
                               withcorr = FALSE,
                               corr = list(),
                               corr.prob = 0.5,
                               tree = list(),
                               de.center = 0,
                               alpha_mean = 0.3, 
                               alpha_sd = 0.002,
                               lenslope=0.02, 
                               nbins=20, 
                               amp_bias_limit = c(-0.2, 0.2),
                               rate_2PCR=0.8, 
                               nPCR1=16, 
                               nPCR2=10, 
                               LinearAmp=F, 
                               LinearAmp_coef=2000, 
                               depth_mean = 100000, 
                               depth_sd = 1500,
                               dirname = "",
                               trials = 1
         ))

#========Generics=========#
#' @rdname setParam
setGeneric("setParam", function(object, name, value) {
  standardGeneric("setParam")
})

#' @rdname setParams
setGeneric("setParams", function(object, update = NULL, ...) {
  standardGeneric("setParams")
})

#' @rdname getParam
setGeneric("getParam", function(object, name) {
  standardGeneric("getParam")
  })

#' @rdname getParams
setGeneric("getParams", function(object, names) {
  standardGeneric("getParams")
})

#' @rdname setParamUnchecked
setGeneric("setParamUnchecked", function(object, name, value) {
  standardGeneric("setParamUnchecked")
})


#' @rdname setParamsUnchecked
setGeneric("setParamsUnchecked", function(object, update = NULL, ...) {
  standardGeneric("setParamsUnchecked")
})


#' @rdname expandParams
setGeneric("expandParams", function(object, vectors, n) {
  standardGeneric("expandParams")
})

#======== Methods=========#
#' Get a parameter
#'
#' Get a parameter value from a Params object.  
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' 
#' @param object Params object to get values from.
#' @param name vector of names of the parameters to get.
#'
#' @return List with the values of the selected parameters.
#' @examples
#' params <- newescoParams()
#' getParam(params, c("nGenes"))
#' @rdname getParam
#' @export
#' @importFrom methods slot
setMethod("getParam", "Params", function(object, name) {
  slot(object, name)
})

#' Get parameters
#'
#' Get multiple parameter values from a Params object. 
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' @param object Params object to get values from.
#' @param names vector of names of the parameters to get.
#'
#' @return List with the values of the selected parameters.
#' @examples
#' params <- newescoParams()
#' getParams(params, c("nGenes", "nCells"))
#' @rdname getParams
#' @export
#' @importFrom methods slot
setMethod("getParams", "Params", function(object, names) {
  checkmate::assertClass(object, classes = "Params")
  checkmate::assertCharacter(names, min.len = 1, any.missing = FALSE)
  params.list <- lapply(names, getParam, object = object)
  names(params.list) <- names
  return(params.list)
}
)

#' Set a parameter
#'
#' Function for setting a parameter value. 
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' 
#' @param object object to set parameter in.
#' @param name name of the parameter to set.
#' @param value value to set the parameter to.
#'
#' @return Object with new parameter value.
#'
#' @examples
#' params <- newescoParams()
#' setParam(params, "nGenes", 100)
#' 
#' @rdname setParam
#' @export
#' @importFrom methods slot<- validObject
setMethod("setParam", "Params", function(object, name, value) {
  checkmate::assertString(name)
  slot(object, name) <- value
  validObject(object)
  return(object)
})

#' Set mutiple parameters
#'
#' Function for setting multiple parameters value.  
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' 
#' @param object Params object to set parameters in.
#' @param update list of parameters to set where \code{names(update)} are the
#'        names of the parameters to set and the items in the list are values.
#' @param ... additional parameters to set. These are combined with any
#'        parameters specified in \code{update}.
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. If the same
#' parameter is specified multiple times it will be set multiple times.
#' Parameters can be specified using a list via \code{update} (useful when
#' collecting parameter values in some way) or individually (useful when setting
#' them manually), see examples.
#'
#' @return Params object with updated values.
#'
#' @examples
#' params <- newescoParams()
#' params
#' # Set individually
#' params <- setParams(params, nGenes = 1000, nCells = 50)
#' params
#' # Set via update list
#' params <- setParams(params, list(mean.rate = 0.2, mean.shape = 0.8))
#' params
#'
#' @rdname setParams
#' @export
setMethod("setParams", "Params", function(object, update = NULL, ...) {
  checkmate::assertClass(object, classes = "Params")
  checkmate::assertList(update, null.ok = TRUE)
  update <- c(update, list(...))
  if (length(update) > 0) {
    for (name in names(update)) {
      value <- update[[name]]
      object <- setParam(object, name, value)
    }
  }
  return(object)
})


#' Set a parameter UNCHECKED
#'
#' Function for setting parameter values without checking for validity. 
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' 
#' @param object object to set parameter in.
#' @param name name of the parameter to set.
#' @param value value to set the parameter to.
#'
#' @return Object with new parameter value.
#' @rdname setParamUnchecked
#' @importFrom methods slot<-
setMethod("setParamUnchecked", "Params", function(object, name, value) {
  checkmate::assertString(name)
  slot(object, name) <- value
  return(object)
})

#' Set parameters UNCHECKED
#'
#' Set multiple parameters in a Params object without checking for validity.
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' 
#' @param object Params object to set parameters in.
#' @param update list of parameters to set where \code{names(update)} are the
#'        names of the parameters to set and the items in the list are values.
#' @param ... additional parameters to set. These are combined with any
#'        parameters specified in \code{update}.
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. If the same
#' parameter is specified multiple times it will be set multiple times.
#' Parameters can be specified using a list via \code{update} (useful when
#' collecting parameter values in some way) or individually (useful when setting
#' them manually), see examples. 
#'
#' @return Params object with updated values.
#' @rdname setParamsUnchecked
setMethod("setParamsUnchecked", "Params", function(object, update = NULL, ...) {
  
  checkmate::assertClass(object, classes = "Params")
  checkmate::assertList(update, null.ok = TRUE)
  
  update <- c(update, list(...))
  
  if (length(update) > 0) {
    for (name in names(update)) {
      value <- update[[name]]
      object <- setParamUnchecked(object, name, value)
    }
  }
  return(object)
})

#' Expand values of parameters 
#'
#' Expand the values of multiple parameters by repititions in a Params object.
#' This function is similar with that in \code{\link[splatter]{splatter}}.
#' 
#' @param object Params object to set parameters in.
#' @param vectors list of parameters to reset values,
#'  where \code{names(update)} are the
#'        names of the parameters to set and the items 
#'        in the list are values.
#' @param n an integer indicates the target number of repititions.
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. 
#' If the same parameter is specified multiple times 
#' it will be set multiple times.
#'
#' @return object with new parameters.
#' @rdname expandParams
setMethod("expandParams", "Params", function(object, vectors, n) {
  update <- list()
  for (parameter in vectors) {
    value <- getParam(object, parameter)
    if (length(value) == 1) {
      update[[parameter]] <- rep(value, n)
    }
  }
  object <- setParamsUnchecked(object, update)
  return(object)
})


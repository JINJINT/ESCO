#' Initiate a new escoParam object
#' 
#' To initiate a new object for storing esco simulation parameters and process.
#' The default parameters are set as in \code{escoParams}.
#' 
#' @param ... any additional parameter settings to override what is provided in
#'        \code{\link{params}}.
#'        
#' @importFrom methods new
#' @return Object with new parameter value.
#' @export
#' @examples 
#' para = newescoParams()
#' @rdname newParams
newescoParams <- function(...) {

    params <- new("escoParams")
    params <- setParams(params, ...)

    return(params)
}

#' Check the validaty of esco parameter object
#' 
#' @import checkmate

setValidity("escoParams", function(object) {
    object <- expandParams(object)
    v <- getParams(object, c(slotNames(object)))

    nBatches <- v$nBatches
    nGroups <- v$nGroups
    checks <- c(nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                mean.rate = checkNumber(v$mean.rate, lower = 0),
                mean.shape = checkNumber(v$mean.shape, lower = 0),
                lib.loc = checkNumber(v$lib.loc),
                lib.scale = checkNumber(v$lib.scale, lower = 0),
                lib.norm = checkFlag(v$lib.norm),
                out.prob = checkNumber(v$out.prob, lower = 0, upper = 1),
                out.facLoc = checkNumber(v$out.facLoc),
                out.facScale = checkNumber(v$out.facScale, lower = 0),
                nGroups = checkInt(v$nGroups, lower = 1),
                group.prob = checkNumeric(v$de.prob, lower = 0, upper = 1,
                                          len = nGroups),
                de.prob = checkNumeric(v$de.prob, lower = 0, upper = 1,
                                       len = nGroups),
                de.downProb = checkNumeric(v$de.downProb, lower = 0, upper = 1,
                                           len = nGroups),
                de.facLoc = checkNumeric(v$de.facLoc, len = nGroups),
                de.facScale = checkNumeric(v$de.facScale, lower = 0,
                                           len = nGroups),
                bcv.common = checkNumber(v$bcv.common, lower = 0),
                bcv.df = checkNumber(v$bcv.df, lower = 0),
                dropout.type = checkCharacter(v$dropout.type,
                                              any.missing = FALSE),
                dropout.mid = checkNumeric(v$dropout.mid, finite = TRUE,
                                           any.missing = FALSE, min.len = 1),
                dropout.shape = checkNumeric(v$dropout.shape, finite = TRUE,
                                             any.missing = FALSE, min.len = 1),
                seed = checkInt(v$seed, lower = 0))


    # Check group.prob sums to 1
    if (sum(v$group.prob) != 1) {
        checks <- c(checks, "group.probs must sum to 1")
    }
     
    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

#' Set an entry esco parameter object
#' @seealso 
#' \code{\link[splatter]{setParam}}
#' @rdname setParam
setMethod("setParam", "escoParams",function(object, name, value) {
    checkmate::assertString(name)
    
    if (name == "group.prob") {
        object <- setParamUnchecked(object, "nGroups", length(value))
    }

    object <- callNextMethod()

    return(object)
})


#' Set multiple entries esco parameter object
#' \code{\link[splatter]{setParams}}
#' @rdname setParams
setMethod("setParams", "escoParams", function(object, update = NULL, ...) {
    
    checkmate::assertClass(object, classes = "escoParams")
    checkmate::assertList(update, null.ok = TRUE)
    
    update <- c(update, list(...))
    
    object <- callNextMethod(object, update)
    
    return(object)
})

#' Expand multiple entries esco parameter object
#' @rdname expandParams
setMethod("expandParams", "escoParams", function(object) {
    
    n <- getParam(object, "nGroups")
    
    vectors <- c("de.prob", "de.downProb", "de.facLoc", "de.facScale")
    
    object <- callNextMethod(object, vectors, n)
    
    return(object)
})



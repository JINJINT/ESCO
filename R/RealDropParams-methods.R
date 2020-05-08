#' @rdname newParams
#' @importFrom methods new
#' @export
newRealDropParams <- function(...) {
  
  params <- new("RealDropParams")
  params <- setParams(params, ...)
  
  return(params)
}

#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
#' checkFlag
setValidity("RealDropParams", function(object) {
  
  object <- expandParams(object)
  v <- getParams(object, c(slotNames(object)))
  
  checks <- c(dropout.mid = checkNumeric(v$dropout.mid, finite = TRUE,
                                         any.missing = FALSE, min.len = 1),
              dropout.shape = checkNumeric(v$dropout.shape, finite = TRUE,
                                           any.missing = FALSE, min.len = 1),
              seed = checkInt(v$seed, lower = 0))
  
  
  if (all(checks == TRUE)) {
    valid <- TRUE
  } else {
    valid <- checks[checks != TRUE]
    valid <- paste(names(valid), valid, sep = ": ")
  }
  
  return(valid)
})

#' @rdname setParam
setMethod("setParam", "RealDropParams", function(object, name, value) {
  checkmate::assertString(name)

  object <- callNextMethod()
  
  return(object)
})

#' @importFrom methods callNextMethod
setMethod("show", "RealDropParams", function(object) {
  
  pp <- list("Dropout:"        = c("(Midpoint)"     = "dropout.mid",
                                   "(Shape)"        = "dropout.shape"))
  
  callNextMethod()
  showPP(object, pp)
})

#' @rdname expandParams
setMethod("expandParams", "RealDropParams", function(object) {

  return(object)
})

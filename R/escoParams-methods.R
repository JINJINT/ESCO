#' @rdname newParams
#' @importFrom methods new
#' @export
newescoParams <- function(...) {

    params <- new("escoParams")
    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric checkCharacter
#' checkFlag
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

    # # Check batchCells matches nCells, nBatches
    # if (v$nCells != sum(v$batchCells) || nBatches != length(v$batchCells)) {
    #     checks <- c(checks,
    #                 "nCells, nBatches and batchesCells are not consistent")
    # }

    # Check group.prob sums to 1
    if (sum(v$group.prob) != 1) {
        checks <- c(checks, "group.probs must sum to 1")
    }

    # # Check path.from
    # if (!(0 %in% v$path.from)) {
    #    checks <- c(checks, path.from = "origin must be specified in path.from")
    # } else if (any(v$path.from == seq_len(nGroups))) {
    #     checks <- c(checks, "path cannot begin at itself")
    # }

    # # Check dropout type
    # if (!(v$dropout.type %in%
    #       c("none", "experiment", "batch", "group", "cell"))) {
    #     checks <- c(checks,
    #                 paste("dropout.type must be one of: 'none', 'experiment',",
    #                       "'batch', 'group', 'cell'"))
    # }
    # 
    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

#' @rdname setParam
setMethod("setParam", "escoParams",function(object, name, value) {
    checkmate::assertString(name)

    # if (name == "nCells" || name == "nBatches") {
    #     stop(name, " cannot be set directly, set batchCells instead")
    # }

    #if (name == "nGroups") {
        
        #stop(name, " cannot be set directly, set group.prob instead")
    #}

    # if (name == "batchCells") {
    #     object <- setParamUnchecked(object, "nCells", sum(value))
    #     object <- setParamUnchecked(object, "nBatches", length(value))
    # }
    # 
    
    if (name == "group.prob") {
        object <- setParamUnchecked(object, "nGroups", length(value))
    }


    object <- callNextMethod()

    return(object)
})

#' @importFrom methods callNextMethod
setMethod("show", "escoParams", function(object) {

    pp <- list("Mean:"           = c("(Rate)"         = "mean.rate",
                                     "(Shape)"        = "mean.shape",
                                     "[Method]"       = "mean.method"),
                                     #"(Density)"      = "mean.dens"),
               "Exprs outliers:" = c("(Probability)"  = "out.prob",
                                     "(Location)"     = "out.facLoc",
                                     "(Scale)"        = "out.facScale"),
               "Groups:"         = c("[Groups]"       = "nGroups",
                                     "(Group Probs)"  = "group.prob"),
               "Tree:"          = c( "[Tree Design]"  = "tree",
                                     "[Tree DE mean]"  = "de.center"),
               "Paths:"          = c("[Path Design]"  = "paths.design",
                                     "[Cell Deisgn]"  = "cells.design",
                                     "[Path mean]"    = "paths.means"),
               "Diff expr:"      = c("(DE Prob in all)"  = "deall.prob",
                                     "(DE Prob in each group)"  = "de.prob",
                                     "(Location)"     = "de.facLoc",
                                     "(Scale)"        = "de.facScale"),
               "Library size:"   = c("(Location)"     = "lib.loc",
                                     "(Scale)"        = "lib.scale",
                                     "(Norm)"         = "lib.norm",
                                     "[Method]"       = "lib.method"),
                                     #"(Density)"      = "lib.dens"),
               "BCV:"            = c("(Common Disp)"  = "bcv.common",
                                     "(DoF)"          = "bcv.df"),
               "Corr:"           = c("[Correlation]"  = "withcorr",
                                     "(Correlation list)"  = "corr",
                                     "(Correlation probablity)"  = "corr.prob"
               ),
               "Dropout:"        = c("[Type]"         = "dropout.type",
                                     "(Midpoint)"     = "dropout.mid",
                                     "(Shape)"        = "dropout.shape",
                                     "[Down Mean]"    = "alpha_mean",
                                     "[Down SD]"      = "alpha_sd",
                                     "[Gene length]"  = "lenslope",
                                     "[Bin numbers]"  = "nbins",
                                     "[Amplification]" = "amp_bias_limit",
                                     "[PCR]" = "rate_2PCR",
                                     "[PCR first round]" = "nPCR1",
                                     "[PCR second round]" = "nPCR2",
                                     "[Linear amplification]" = "LinearAmp",
                                     "[Linear amplification coef]" = "LinearAmp_coef",
                                     "[Depth mean]" = "depth_mean",
                                     "[Depth sd]" = "depth_sd"
                                     )
        )
    callNextMethod()
    showPP(object, pp)
})

#' @rdname expandParams
setMethod("expandParams", "escoParams", function(object) {

    n <- getParam(object, "nGroups")

    vectors <- c("de.prob", "de.downProb", "de.facLoc", "de.facScale")

    object <- callNextMethod(object, vectors, n)

    return(object)
})


#' @rdname setParams
setMethod("setParams", "escoParams", function(object, update = NULL, ...) {
    
    checkmate::assertClass(object, classes = "escoParams")
    checkmate::assertList(update, null.ok = TRUE)
    
    update <- c(update, list(...))
    
    update <- bringItemsForward(update, c("nCells", "group.prob"))
    
    object <- callNextMethod(object, update)
    
    return(object)
})

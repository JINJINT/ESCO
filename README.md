# ESCO

Welcome to ``ESCO``! ``ESCO`` is an R package for the simulation of
single-cell RNA sequencing data with special consideration of gene co-expression, built ultilizing the infrastructure in the [`splatter`](https://github.com/Oshlack/splatter) package. 

## Install from Github
This package can be installed with R package devtools:
```{r}
library("devtools")
devtools::install_github("JINJINT/ESCO")
```

## Quick start:

For a simple example to simulate 100 genes and 50 cells of one cell group with gene co-expression:
```{r}
library(ESCO)

#===== start simulation ======#
sim <- escoSimulateSingle(nGenes = 100, nCells = 50, 
                          withcorr = TRUE,
                          verbose = FALSE)

#===== access the data ======#
datalist = list("simulated truth"=assays(sim)$TrueCounts,
                "zero-inflated" = assays(sim)$counts, 
                "down-sampled" = assays(sim)$observedcounts)

#====== plot the data ======#
heatdata(datalist, norm = FALSE, size = 2, ncol = 3)

#====== plot the Gene correlation ======#
# object that saved all simulation configurations
simparams = metadata(sim)$Params 

# object that particularly saved the correlation structure
rholist = slot(simparams,"corr") 

# arrange the true correlation and simulated correlation
corrgenes = rownames(rholist[[1]])
gcnlist = lapply(datalist, function(data)gcn(data, genes = corrgenes))
gcnlist = append(gcnlist, list("given truth" = rholist[[1]]), 0)
heatgcn(gcnlist, size = 3, ncol = 4)
```

For more complicated examples of simulating multiple cell groups and even trees and trajectories with gene co-expression, please check out the [vignettes](https://www.dropbox.com/s/ly1x20c7bommsvi/esco.html?dl=0), which can also be built locally if installed by 
```{r}
devtools::install_github("JINJINT/ESCO", build_vignettes=TRUE)
```

## Reference:
Check out our paper for ESCO here:
[ESCO: single cell expression simulation incorporating gene co-expression.](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab116/6149079?guestAccessKey=64c91aa4-1d5e-42da-92df-678b1b08af79)

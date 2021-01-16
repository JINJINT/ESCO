# ESCO

Welcome to ``ESCO``! ``ESCO`` is an R package for the simulation of
single-cell RNA sequencing data with special consideration of gene co-expression, built ultilizing the infrastructure the [`splatter`](https://github.com/Oshlack/splatter) package. This package formerly known as `SplatterESCO`.

## Install from Github
This package can be installed with R package devtools. First, pull the package with git clone to your working directory. Make sure that you have installed the packages listed in the DESCRIPTION file. 

To install ESCO, run:
```{r}
library("devtools")
devtools::install_github("JINJINT/ESCO")
```

## Quick start:

For a simple example to simulate 100 genes and 50 cells of one cell group with gene co-expression:
```{r}
# start simulation
sim <- escoSimulateSingle(nGenes = 100, nCells = 50, 
                          withcorr = TRUE,
                          verbose = FALSE)

#===== access the data
datalist = list("simulated truth"=assays(sim)$TrueCounts,
                "zero-inflated" = assays(sim)$counts, 
                "down-sampled" = assays(sim)$observedcounts)

#====== plotting the data
heatdata(datalist, norm = FALSE, size = 2, ncol = 3)

#====== plotting the Gene correlation
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
[ESCO: single cell expression simulation incorporating gene co-expression.](https://www.biorxiv.org/content/10.1101/2020.10.20.347211v1)

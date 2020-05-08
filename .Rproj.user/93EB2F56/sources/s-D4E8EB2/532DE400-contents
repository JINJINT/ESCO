
# evaluate gcn
# return error summary
#' @import ggplot2
#' @import utils
#' @import stats
#' @import AUC
#' @import clues
#' @export
evaluation <- function(datalist, genes, cells,  plotdir = NULL,
                       corrname = c("bayescor", "pearson", "spearman", "cosine", "kendall"),
                       CPM2 = FALSE,   
                       errname = c("MSE", "Spectral", "Inf", "L1", "ARI", "Jaccard", "AUC"), 
                       clustk = c(2:14), clustmethod = "hclust"){
    
  if("ARI" %in% errname | "Jaccard" %in% errname){

      if(max(clustk)>length(genes)){
        clustk = length(genes)-1
      }
      if(min(clustk)==1){
        cat("Sorry, requires at least three genes to compute the ARI/ Jaccard/ AUC error...\n")
        return(NULL)
      }
      if("ARI" %in% errname){
        errname = errname[-c(which(errname=="ARI"))]
        errname = c(errname, paste0("ARI_k",clustk))
      }
     if("Jaccard" %in% errname){
       errname = errname[-c(which(errname=="Jaccard"))]
       errname = c(errname, paste0("Jaccard_k",clustk))
     }
  }
    
    gcnerror = array(0, dim = c(length(errname), length(corrname), length(datalist)-1))
    dimnames(gcnerror) = list(errname, corrname, names(datalist)[2:length(datalist)])
    
    for(corr in corrname){
      gcnlist = lapply(datalist, function(data){
        if(sum(is.na(data))==nrow(data)*ncol(data))return(matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes)))
        subdata = data[,cells]
        subgcn = gcn(subdata, genes, CPM  = TRUE, CPM2 = CPM2, name = corr)
        return(subgcn)
      })
      names(gcnlist)  = names(datalist)
      if(corr == "pearson" & !(is.null(plotdir)))heatgcn(gcnlist, dirname = plotdir)
      
      gcn_true = gcnlist[[1]]
      gcn_true[is.na(gcn_true)] = 0
      
      if("AUC" %in% errname){
        truedata = datalist[[1]]
        truedata = 10^6*t(t(truedata)/colSums(truedata))
        data = truedata[genes, cells]
        if(CPM2)data = t(t(data)/colSums(data))
        data[is.na(data)]=0
        
        # parallel
        numCores <- detectCores() -1
        cl <- makeCluster(numCores)
        registerDoSNOW(cl)
        
        combgenes = t(combn(genes, 2))
        total <- nrow(combgenes)
        pb <- progress_bar$new(
          format = "progress = :letter [:bar] :elapsed | eta: :eta", total = total, width = 60)
        progress <- function(n){
          pb$tick(tokens = list(letter = rep("", total)[n]))
        } 
        opts <- list(progress = progress)
        corpval <- foreach(i = 1:nrow(combgenes), .combine = cbind, .options.snow = opts ) %dopar% {
                        ans = stats::cor.test(data[combgenes[i,1], ], data[combgenes[i,2], ], alternative = "two.sided", exact = TRUE)$p.value
                        return(ans)
        }
        stopCluster(cl) 
        
        thred = min(quantile(corpval, 0.3, na.rm = TRUE), 0.05)
        edges = 1*(corpval <= thred)
        edges[is.na(edges)] = 0
        noedges = 1-edges
        edges = as.factor(edges)
        noedges = as.factor(noedges)
      }
      
      for(err in errname){
        compute_err <-function(gcnimp){
           gcnimp[is.na(gcnimp)] = 0
           gcnimp_diff = gcn_true-gcnimp
           gcnimp_diff[is.na(gcnimp_diff)]=0
          
           if(err=="MSE")gcnimp_error = (norm(gcnimp_diff, type = "F")/nrow(gcnimp_diff))^2
           if(err=="Spectral")gcnimp_error = norm(gcnimp_diff, type = "2")
           if(err=="Inf")gcnimp_error = norm(gcnimp_diff, type = "I")
           if(err=="L1")gcnimp_error = norm(gcnimp_diff, type = "O")
           if(err== "AUC"){
            impvec = apply(t(combn(genes, 2)), 1, function(x)abs(gcnimp[x[1], x[2]])) 
            improc = AUC::roc(impvec, edges)
            gcnimp_error = AUC::auc(improc)
           }
           else{
               clk <- gregexpr("[0-9]+", err)
               clk <- as.integer(unlist(regmatches(err, clk)))
               
               errhere <- gregexpr("[A-Za-z]+", err)
               errhere <- unlist(regmatches(err, errhere))[1]
    
               if(clustmethod == "hclust"){
                 d_true = stats::as.dist((1-gcn_true)/2)
                 h_true = stats::hclust(d_true)
                 memb_true <- cutree(h_true, k = clk)
                 d = stats::as.dist((1-gcnimp)/2)
                 h = stats::hclust(d)
                 memb <- cutree(h, k = clk)
               }
               if(clustmethod == "kmeans"){
                 memb_true <- kmeans(gcn_true, clk)$clust
                 memb <- kmeans(gcnimp, clk)$clust
               }
               true_clust = memb_true
               clust = memb
               if(errhere=="ARI")gcnimp_error =  adjustedRand(clust, true_clust)[2]
               if(errhere=="Jaccard")gcnimp_error =  adjustedRand(clust, true_clust)[5]
           }
           return(gcnimp_error)
        }
        gcnerror[which(errname==err), which(corrname==corr), ] = sapply(gcnlist[-c(1)], compute_err)
      }
    }
    return(gcnerror)
  }




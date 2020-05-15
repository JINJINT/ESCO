# compute gcn 
# return gcn matrix
#' @import ggplot2
#' @export
gcn<-function(count, genes = NULL, CPM = TRUE, CPM2 = FALSE, 
              name = c("pearson"), fast = FALSE){
  if(is.null(genes))genes = rownames(count)
  if(length(name)>1)name = "pearson"
  if(name=="bigscale")gcn = as.numeric(compute.network(count[genes,], rownames(count[genes,]), clustering = "recursive", quantile.p = 0.999)$correlations)
  if(CPM){
      count = 10^6*t(t(count)/colSums(count))
  }
  count = count[as.character(genes),]
  if(CPM2)count = 10^6*t(t(count)/colSums(count))
  count[is.na(count)] = 0
  #if(name=="bayescor")gcn = Bayes_Corr_Prior3(count)
  if(name=="pearson")gcn = cor(t(log2(count+1)))
  if(name=="spearman")gcn = cor(t(log2(count+1)), method = "spearman")
  if(name=="cosine")gcn = cosine(t(count))
  #if(name=="copula")gcn = estimate_copula(log2(count+1), fast)
  if(name=="kendall")gcn = cor(t(log2(count+1)), method = "kendall")
  gcn[is.na(gcn)]=0
  return(gcn)
}

# compute the erros between two gcn
# return the error value
gcn_error<-function(gcnlist, 
                    method = c("MSE", "Spectral", "Inf", "L1", "GO-score", "ARI", "Jaccard"),
                    clustk = 9, clustmethod = "hclust"){
  gcn_true = gcnlist[[1]]
  gcn_true[is.na(gcn_true)] = 0
  gcnerror = c()
  for(i in 2:length(gcnlist)){
    gcn = gcnlist[[i]]
    gcn[is.na(gcn)] = 0
    gcn_diff = gcn_true-gcn
    gcn_diff[is.na(gcn_diff)]=0
    
    if(method=="MSE")gcn_error = (norm(gcn_diff, type = "F")/nrow(gcn_diff))^2
    if(method=="Spectral")gcn_error = norm(gcn_diff, type = "2")
    if(method=="Inf")gcn_error = norm(gcn_diff, type = "I")
    if(method=="L1")gcn_error = norm(gcn_diff, type = "O")
    if(method=="GO-score"){
      load("/Users/jinjin/ADAproject/CoExpNets/bin/run_GBA.Rdata")
      load("/Users/jinjin/ADAproject/CoExpNets/data/GO.human.Rdata")
      source("/Users/jinjin/ADAproject/CoExpNets/bin/helper_functions.r")
      gcn_score = run_GBA(gcn, GO.labels)
      gcn_true_score = run_GBA(gcn_true, GO.labels)
      gcn_error = gcn_true_score - gcn_score
    }
    if(method %in% c("ARI", "Jaccard")){
      if(clustmethod == "hclust"){
        d_true = stats::as.dist((1-gcn_true)/2)
        h_true = stats::hclust(d_true)
        memb_true <- cutree(h_true, k = clustk)
        d = stats::as.dist((1-gcn)/2)
        h = stats::hclust(d)
        memb <- cutree(h, k = clustk)
      }
      if(clustmethod == "kmeans"){
        memb_true <- kmeans(gcn_true, clustk)$clust
        memb <- kmeans(gcn, clustk)$clust
      }
      true_clust = memb_true
      clust = memb
      if(method=="ARI")gcn_error =  adjustedRand(clust, true_clust)[2]
      if(method=="Jaccard")gcn_error =  adjustedRand(clust, true_clust)[5]
    }
    gcnerror[i-1] = gcn_error
  }
  return(gcnerror)
}




# normsc<-function(count){
#   DataNorm = SCnorm(Data = count, Conditions = rep(c(1), each= ncol(count)), FilterCellNum = 10, NCores=6, useZerosToScale=TRUE)
#   ans = results(DataNorm)
#   return(ans)
# }
# 
# 
# ecdf_trans<-function(vec){
#   f = ecdf(vec)
#   return(sapply(vec, function(v)f(v)))
# }
# 
# 
# modi_ecdf<-function(vec, delta){
#   v = ecdf(vec)(vec)
#   vv = (1-delta)*(v>1-delta) + v*(v<=1-delta && v>=delta) + delta*(v<delta)
#   return(vv)
# }
# 
# estimate_copula<-function(data, fast = FALSE ) { 
#   corr_mat <- matrix(1, nrow= nrow(data), ncol = nrow(data))
#   data_col <- matrix(c(0,0), nrow=ncol(data), ncol = 2)
#   if(fast){
#     trans_data = apply(data, 1, ecdf_trans)
#     corr_mat = cor(trans_data)
#   }
#   else{
#     for (i in 1:(nrow(data)-1)){
#       progress(i, progress.bar = TRUE)
#       Sys.sleep(0.01)
#       if (i == (nrow(data)-1)) cat("Done!\n")
#       for (j in (i+1):nrow(data)){
#         data_col[,1]=t(data[i,])
#         data_col[,2]=t(data[j,])
#         m = pobs(data_col)
#         
#         cop = normalCopula(dim = 2)
#         fit = fitCopula(cop, m, method = 'ml')
#         corr_mat[i,j] = coef(fit)[1]
#         corr_mat[j,i] = corr_mat[i,j] 
#       }
#     }
#   }
#   
#   return(corr_mat)
# }





#' Compute the gene co-expression network 
#' 
#' compute the gene co-expression netowrk using different correlation metric
#' @param count a gene by cell matrix containing the expression count.
#' @param genes a vector of names of selected genes, default is NULL, indicating choosing all the genes.
#' @param CPM whether to use Counts Per Millon normalization on the whole data sets or not.
#' @param CPM2 whether to use Counts Per Millon normalization on the selected genes or not.
#' @param name a string indicating what type of correlation metric ("pearson", "spearman", "kendall","cosine") to use, 
#'        default options is "pearson". 
#' @return a correlation matrix
#' @rdname escogcn
#' @examples 
#' data = matrix(rnorm(100),20,5)
#' gcndata = gcn(data)
#' @export
gcn<-function(count, genes = NULL, CPM = TRUE, CPM2 = FALSE, 
              name = "pearson"){
  if(is.null(genes))genes = rownames(count)
  if(CPM){
      count = 10^6*t(t(count)/colSums(count))
  }
  count = count[as.character(genes),]
  if(CPM2)count = 10^6*t(t(count)/colSums(count))
  count[is.na(count)] = 0
  if(name=="pearson")gcn = cor(t(log2(count+1)))
  if(name=="spearman")gcn = cor(t(log2(count+1)), method = "spearman")
  if(name=="cosine")gcn = cosine(log2(count+1))
  if(name=="kendall")gcn = cor(t(log2(count+1)), method = "kendall")
  gcn[is.na(gcn)]=0
  return(gcn)
}

cosine<-function(DF){
  Matrix <- as.matrix(DF)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  return(sim)
}

#' Compute the estiamtion error of gene co-expression networks
#' 
#' Compute the error between the true and estimated gene co-expression network using different error measure
#' @param gcnlist a list of gene coexpression networks of the same set of genes, 
#'        where the first is taken as the truth, and the rest are all estimations.
#' @param method a vector of error measure names, possible choices are mean square error "MSE", 
#'        Spectral norm "Spectral", infinite norm "Inf", l1 norm, "L1", Adjusted rand index \code{\link{adjustedRand}} with 
#'        \code{clusk} numbers of clusters "ARI", Jaccard index with 
#'        \code{clusk} numbers of clusters "Jaccard".
#' @param clustk the number of clusters used in "ARI" and "Jaccard" error measure. Default is 9.
#' @param clustmethod what clustering method to use. Possible choices are hirarchical clustering "hclust" \code{\link{hclust}}, 
#'        and "kmeans" \code{\link{kmeans}}. Default is "hclust".
#' @return a matrix of error levels, where each row indicates the given different estimation of gcn, 
#'         and each column indicates different error measure.
#' @rdname escogcnerr
#' @importFrom stats as.dist hclust
#' @importFrom clues adjustedRand
#' @examples 
#' data1 = matrix(rnorm(100),20,5)
#' data2 = matrix(rnorm(100),20,5)
#' gcnlist = list(cor(data1),cor(data2))
#' gcnerr = gcn_error(gcnlist)
#' @export
gcn_error<-function(gcnlist, 
                    method = c("MSE", "Spectral", "Inf", "L1", "ARI", "Jaccard"),
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
    if(method %in% c("ARI", "Jaccard")){
      if(clustmethod == "hclust"){
        d_true = as.dist((1-gcn_true)/2)
        h_true = hclust(d_true)
        memb_true <- cutree(h_true, k = clustk)
        d = as.dist((1-gcn)/2)
        h = hclust(d)
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







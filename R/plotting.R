#' Visualize a set of gene expression matrix
#'
#' Heatmap plot a set of gene expression matrix of the same set of gene and cells with proper annotations
#'
#' @param datalist a list of gene expression matrixes of size p by n, (where rows are of the same set of genes and columns are of the same set of cells).
#' @param dirname a string of directory names to save the plots, default as NULL, that is not saving but directly showing the plots.
#' @param genes a vector contains names of the genes to be plotted.
#' @param cellinfo a dataframe contains a column named 'newcelltype' of length n, where each entries corresponds to 
#'        the cell group identity that a cell belongs to. This information will be used for annotating the final heatmap.
#'        Default is NULL, that is no cell annotation.
#' @param geneinfo a dataframe contains a column named 'newcelltype' of length p, where each entries corresponds to 
#'        the cell group identity that a gene marks. This information will be used for annotating the final heatmap.
#'        Default is NULL, that is no gene annotation.
#' @param rowv a vector contains the ordering of the genes for all data matrix in the final heatmap.
#'        Default is FALSE, that the ordering resulted from the clustering (\code{hclust}) result of the first data matrix in datalist.
#' @param colv a vector contains the ordering of the cells for all data matrix in the final heatmap.
#'        Default is FALSE, that the ordering resulted from the clustering (\code{hclust}) result of the first data matrix in datalist.
#' @param maxdata the maximun cutoff of the correlation, default is NULL, that is no cutoff.
#' @param mindata the minmum cutoff of the correlation, default is 0.
#' @param color what set of color panel to use, default is "YlGnBu:100".
#' @param extrainfo a string of extra information to saved in the final plot filename, default as "".
#' @param ncol the number of columns in the combined plots, default as 3.
#' @param size the size of the cellwidth and cellheight in heatmap, default as 3.
#' @param width the width of the final pdf plot, default as 8.
#' @param height the height of the final pdf plot, default as 8.
#' @param log whether take the log transform of the data
#' @param norm whether take the counts per million normalization of the data
#' @importFrom NMF aheatmap
#' @import viridis
#' @import RColorBrewer
#' @examples 
#' data = matrix(rnorm(100),20,5)
#' # heatdata(list(data))
#' @rdname escoheatdata
#' @export
heatdata<-function(datalist, dirname = NULL, genes = NULL, cellinfo = NULL, 
                   rowv = FALSE, colv = FALSE, geneinfo = NULL,
                   log = TRUE, norm = TRUE, 
                   maxdata = NULL, mindata = 0, ncol = 3, size = 3, width = 8, height = 8, 
                   color = "YlGnBu:100", extrainfo = ""){

  typecell = NA
  typegene = NA
  cellcolors = NA
  genecolors = NA

  if(norm){
    for(i in 1:length(datalist)){
      data = datalist[[i]]
      datalist[[i]] = 10^6*t(t(data)/colSums(data))
    }
    norm = "norm"
  }
  else norm = ""

  if(is.null(genes))genes = as.character(rownames(datalist[[1]]))
  genegroups = c()
  cellgroups = c()
  if(!is.null(cellinfo))cellgroups= levels(cellinfo$newcelltype)
  if(!is.null(geneinfo))genegroups= levels(geneinfo$newcelltype)
  colors = brewer.pal(n = length(unique(c(genegroups, cellgroups))), name = "Set1")

  if(!is.null(cellinfo)){
    typecell = data.frame("Cell Group" = cellinfo$newcelltype)
    cellcolors =  colors[1:length(levels(cellinfo$newcelltype))]
    colors = colors[-(1:length(levels(cellinfo$newcelltype)))]
  }

  if(!is.null(geneinfo)){
    typegene = data.frame("Gene Type" = geneinfo$newcelltype)
    genecolors = rep("", length(levels(geneinfo$newcelltype)))
    genecolors[match(levels(cellinfo$newcelltype), levels(geneinfo$newcelltype))] = cellcolors
    genecolors[which(genecolors=="")] = colors[1:length(genecolors[which(genecolors=="")])]
    genecolors[which(levels(geneinfo$newcelltype)=="None")]="gray70"
  }

  colr = list("Cell.Group" = cellcolors, "Gene.Type" = genecolors)
  
  if(!is.null(dirname)){

     a = aheatmap(log2(datalist[[1]][genes,]+1), Rowv = rowv, Colv = colv,
           color = color, breaks = seq(0, maxfunc(maxdata, log2(datalist[[1]][genes,]+1)), length.out = 101),
           annRow = typegene, annCol = typecell, annColors = colr, width = width, height = height, filename = paste0(dirname, norm, "data", names(datalist)[1],extrainfo, ".png"))
  
     rowv = a$rowInd
     colv = a$colInd

     datalist = datalist[-1]
     if(length(datalist)>0){
        for(i in 1:(length(datalist))){
          aheatmap(log2(datalist[[i]][genes,] + 1),  Rowv = rowv, Colv = colv, revC = FALSE,
               color = color, breaks = seq(mindata, maxfunc(maxdata,log2(datalist[[i]][genes,]+1)), length.out = 101),
               annRow = typegene, annCol = typecell, annColors = colr, width = width, height = height,
               filename = paste0(dirname, norm, "data", names(datalist)[i],extrainfo,".png"))
      }
    }
  }
  
  else{
    par(mfrow=c(ceiling(length(datalist)/ncol), ncol))
    legend = FALSE
    if(length(datalist)==1)legend = TRUE
    a = aheatmap(log2(datalist[[1]][genes,]+1), Rowv = rowv, Colv = colv,
                 color = color, breaks = seq(0, maxfunc(maxdata, log2(datalist[[1]][genes,]+1)), length.out = 101),
                 annRow = typegene, annCol = typecell, annColors = colr, legend = legend, annLegend = legend, cellwidth = size, cellheight = size,  main = names(datalist)[1])
    
    rowv = a$rowInd
    colv = a$colInd
    if(is.null(maxdata))maxdata = max(log2(datalist[[1]][genes,]+1))
    
    datalist = datalist[-1]
    if(length(datalist)>0){
      for(i in 1:(length(datalist))){
        legend = FALSE
        if(i==length(datalist))legend = TRUE
        
        aheatmap(log2(datalist[[i]][genes,] + 1),  Rowv = rowv, Colv = colv, revC = FALSE,
                 color = color, breaks = seq(mindata, maxfunc(maxdata,log2(datalist[[i]][genes,]+1)), length.out = 101),
                 annRow = typegene, annCol = typecell, annColors = colr, legend = legend, annLegend = legend, cellwidth = size, cellheight = size, main = names(datalist)[i])
      }
    }
    
  }

}

#' Visualize a set of gene correlation matrix
#'
#' Heatmap plot a set of gene correlation matrix of the same set of genes with proper annotations
#' 
#' @param gcnlist a list of gene correlation matrixes of size p by p, where the rows and columns are of the same set of genes.
#' @param dirname a string of directory names to save the plots, default as NULL, that is not saving but directly showing the plots.
#' @param geneinfo a dataframe contains a column named 'newcelltype' of length p, where each entries corresponds to 
#'        the cell group identity that a gene marks. This information will be used for annotating the final heatmap.
#'        Default is NULL, that is no annotation.
#' @param ord a vector contains the ordering of the genes for all correlation matrix in the final heatmap.
#'        Default is NULL, that the ordering resulted from the clustering (\code{hclust}) result of the first gene correlation matrix in gcnlist.
#' @param maxgcn the maximun cutoff of the correlation, default is 1.
#' @param mingcn the minmum cutoff of the correlation, default is -1.
#' @param abs whether take the absulute value of the correlation matrix or not.
#' @param color what set of color panel to use, default is "-RdBu:100".
#' @param extrainfo a string of extra information to saved in the final plot filename, default as "".
#' @param ncol the number of columns in the combined plots, default as 4.
#' @param size the size of the cellwidth and cellheight in heatmap, default as 2.
#' @return heatmap plots showing immediately or pdf files saved to desinated directory.
#' @importFrom NMF aheatmap
#' @importFrom stats as.dist hclust
#' @import viridis
#' @export
#' @examples 
#' data1 = matrix(rnorm(100),20,5)
#' data2 = matrix(rnorm(100),20,5)
#' gcnlist = lapply(list(data1,data2),function(data)gcn(data))
#' # heatgcn(gcnlist)
#' @rdname escoheatgcn
heatgcn<-function(gcnlist, dirname = NULL, geneinfo = NULL,
                  ord = NULL,
                  maxgcn = 1, mingcn = -1, abs = FALSE,
                  color="-RdBu:100", extrainfo = NULL, ncol = 4, size = 2){

  if(sum(mean(sapply(gcnlist, nrow))!=sapply(gcnlist, nrow))>0){
    for(i in 2:length(gcnlist)){
      gcnold = gcnlist[[i]]
      gcnlist[[i]] = matrix(0, nrow(gcnlist[[1]]), ncol(gcnlist[[1]]))
      dimnames(gcnlist[[i]]) = dimnames(gcnlist[[1]])
      rows = rownames(gcnold)[which(rownames(gcnold)%in%rownames(gcnlist[[1]]))]
      gcnlist[[i]][rows,  rows] = gcnold[rows, rows]
    }
  }

  if(is.null(ord)){
    d = as.dist((1-gcnlist[[1]])/2)
    h = hclust(d)
    ord = h$order
    #memb <- cutree(h, k = clustk)
    if(!is.null(geneinfo)){
      #type = data.frame(TypeC = as.factor(memb[ord]), TypeG = geneinfo$newcelltype[ord])
      type = data.frame("Gene Type" = geneinfo$newcelltype[ord])
      genecolors =  brewer.pal(n = length(levels(as.factor(geneinfo$newcelltype))), name = "Set1")
      genecolors[which(levels(geneinfo$newcelltype)=="None")]="gray70"
      #clustcolors = viridis_pal(option = "D")(clustk)
      #colr = list(TypeC =  clustcolors, TypeG = genecolors)
      colr = list("Gene.Type" = genecolors)
    }
    else{
      type = NA
      colr  = NA
      }
  }
  else{
    type = data.frame("Gene Type" = geneinfo$newcelltype[ord])
    genecolors =  brewer.pal(n = length(levels(as.factor(geneinfo$newcelltype))), name = "Set1")
    genecolors[which(levels(geneinfo$newcelltype)=="None")]="gray70"
    colr = list("Gene.Type" = genecolors)
  }

  if(abs&!is.null(mingcn))mingcn = 0
  if(!is.null(dirname)){
    for(i in 1:length(gcnlist)){
      gcn = gcnlist[[i]]
    
      aheatmap(gcn[ord, ord],  Rowv = NA, Colv = NA,
               color = color, breaks = seq(mingcn, maxgcn,length.out = 101),
               annCol = type, annColors = colr, cellwidth = size, cellheight = size, 
               filename = paste0(dirname, "gcn_", names(gcnlist)[i], extrainfo, ".pdf") )
    }
  }   
  else{
    par(mfrow = c(ceiling(length(gcnlist)/ncol), ncol))
    for(i in 1:length(gcnlist)){
      gcn = gcnlist[[i]]
      legend = FALSE
      if(i==length(gcnlist))legend = TRUE
      aheatmap(gcn[ord, ord],  Rowv = NA, Colv = NA,
               color = color, breaks = seq(mingcn, maxgcn, length.out = 101),
               annCol = type, cellwidth = size, cellheight = size, legend = legend, annLegend = legend, annColors = colr, main = paste0(names(gcnlist)[i], " ", extrainfo))
      
    }
  }
}


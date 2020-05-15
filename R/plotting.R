
#' @export
maxfunc <- function(x, data){
  if(is.null(x))return(max(data))
  else return(x)
}

#' @export
minfunc <- function(x, data){
  if(is.null(x))return(min(data))
  else return(x)
}



#' @import NMF
#' @import viridis
#' @import RColorBrewer
#' @export
# plot the gene expressions
heatdata<-function(datalist, dirname = NULL, genes = NULL, cellinfo = NULL, rowv = FALSE, colv = FALSE, geneinfo = NULL,
                   log = TRUE, norm = TRUE, order = NULL,
                   maxdata = NULL, mindata = 0, ncol = 3, size = 3, width = 8, height = 8, color = "YlGnBu:100", extrainfo = ""){

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


#' @import NMF
#' @import viridis
#' @export
# plot the gene correlation matrix
heatgcn<-function(gcnlist, dirname = NULL, geneinfo = NULL,
                  ord = NULL,
                  maxgcn = 1, mingcn = -1, abs = FALSE,
                  color="-RdBu:100", extrainfo = NULL, breaks = 0, ncol = 4, size = 2){

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
    d = stats::as.dist((1-gcnlist[[1]])/2)
    h = stats::hclust(d)
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
    # else{
    #   type = data.frame(TypeC = as.factor(memb[ord]))
    #   clustcolors = viridis_pal(option = "D")(clustk)
    #   colr = list(TypeC =  clustcolors)
    # }
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


#' @import NMF
#' @import viridis
#' @import ggplot2
#' @export
# compare the gene correlation matrix
gcn_compare<-function(gcnavg, gcnstd, xlab="dropout level", normname, corname, colors, main, dirname){
  mdata = melt(gcnavg)
  mdatavar = melt(gcnstd)
  mdata$sd = mdatavar$value
  mdata$Var1 = as.numeric(mdata$Var1)
  mdata$Var3 = as.factor(mdata$Var2%ni%c("Raw", "EnImpute"))
  mdata$Var4 = as.factor(0 + 1*(mdata$Var2=="Raw") + 2*(mdata$Var2=="EnImpute"))
  width = 0.01*(mdata$Var1[2]-mdata$Var1[1])
  ggplot(data = mdata, aes(x = Var1, y = value, color = as.factor(Var2))) +
        geom_line(aes(alpha = Var3, size = Var3, linetype = Var4))+
        geom_point(aes(alpha = Var3, size = Var3)) +
        scale_alpha_discrete(range = c(1, 1))+
        scale_size_manual(values = c(1, 0.4))+
        scale_linetype_manual(values=c("solid", "solid", "twodash")) +
        geom_errorbar(data = mdata, aes(ymin = value - sd, ymax = value + sd, color = Var2), width = width)+
        #scale_shape_manual(values=c(10, 3)) +
        labs(x = xlab, y = normname, color = "Imputation \n methods", title =paste(main, corname)) +
        scale_color_manual(labels = colnames(gcnavg), values=colors) +
        theme_bw() +
        theme(axis.text.y = element_text(size=13),
              axis.text.x = element_text(size=13),
              axis.title.x = element_text(size=13),
              axis.title.y = element_text(size=13),
              plot.title=element_text(hjust=0,size=16)
        )+
        scale_x_continuous(breaks = mdata$Var1)
  ggsave(filename = paste0(dirname, main, corname, ".png"), width = 6, height = 4, units = "in")
}

#' @export
data_roc<-function(datalist, subgenes, xlim){
  truedata = datalist[[1]]
  data = t(t(truedata)/colSums(truedata))
  data = t(t(data)/colSums(data))
  data[is.na(data)] = 0
  corpval = apply(t(combn(subgenes, 2)), 1,
                  function(x)cor.test(data[x[1],], data[x[2],],
                                      alternative = "two.sided", exact = TRUE)$p.value)
  truegcn = gcn(truedata, subgenes, norm = FALSE, name = "pearson")
  truevec = apply(t(combn(subgenes, 2)), 1, function(x)abs(truegcn[x[1], x[2]]))
  edges = 1*(corpval < 0.05)
  print("The proportion of edges is", mean(edges))
  edges = as.factor(edges)
  rocfdr = c()
  rocpow = c()
  imoauc = c()
  for(impdata in datalist[-1]){
    impgcn = gcn(impdata, subgenes, name = "pearson")
    impvec = apply(t(combn(subgenes, 2)), 1, function(x)abs(impgcn[x[1], x[2]]))
    plot(truevec, impvec, main = paste0(dm, "Bag", bag))
    improc = roc(impvec, edges)
    impauc = cbind(auc(improc), 0, 0.2)
    impfdr = improc$fpr
    imppow = improc$tpr
    rocfdr = cbind(rocfdr, impfdr)
    rocpow = cbind(rocpow, imppow)
  }
  colnames(rocfdr) = names(datalist)[-1]
  colnames(rocpow) = names(datalist)[-1]
  gcn_roc(rocfdr, rocpow, xlim = xlim, title = title, colors = c("black", brewer.pal(n = length(datalist)-2, name = "Set1")))
}



#' @export
gcn_roc<-function(fdrvec, powvec, xlab, ylab, title, colors, dirname, xlim = 0.15){
  xdata = melt(fdrvec)
  ydata = melt(powvec)
  data = data.frame(fdr = xdata$value, pow = ydata$value, method = xdata$Var2)
  data$Var3 = as.factor(grepl("Raw", data$method) | grepl("EnImpute", data$method) )
  data$Var4 = as.factor(0 + 1*grepl("Raw", data$method) + 2*grepl("EnImpute", data$method))

  ggplot(data = data, aes(x = fdr, y = pow, color = as.factor(method))) +
          #geom_smooth(aes(size = Var3, linetype = Var4), se = FALSE, span = 10, method = "glm") +
          geom_line(aes(size = Var3, linetype = Var4, alpha = Var3)) +
          scale_alpha_discrete(range = c(1, 1))+
          scale_size_manual(values = c(0.4, 1))+
          scale_linetype_manual(values=c("solid", "solid", "twodash")) +
          #scale_shape_manual(values=c(19, 19, 4)) +
          labs(x = xlab, y = ylab, color = "Imputation \n methods", title =title) +
          scale_color_manual(labels = colnames(fdrvec), values=colors) +
          theme_bw() +
          theme(axis.text.y = element_text(size=13),
                axis.text.x = element_text(size=13),
                axis.title.x = element_text(size=16),
                axis.title.y = element_text(size=16),
                plot.title=element_text(hjust=0,size=16)
          )+
          guides(colour = guide_legend(override.aes = list(size=5)))+
          xlim(0,xlim)+
          ylim(0,1)
  ggsave(filename = paste0(dirname, title,".png"), width = 6, height = 4, units = "in")
}

#' @export
# find gene pairs that most different/similar in two gcn
gene_pairs<- function(gcn1, gcn2 = NULL, low, high){
  if(is.null(gcn2))gcn2 = matrix(0, nrow = nrow(gcn1), ncol = ncol(gcn1))
  diffnet = abs(gcn1 - gcn2)
  genenames = colnames(diffnet)

  diffpairs = which(diffnet>=high, arr.ind = TRUE)
  diffpairs  = subset(diffpairs, diffpairs[,1]!=diffpairs[,2])
  diffgenes1 = genenames[diffpairs[,1]]
  diffgenes2 = genenames[diffpairs[,2]]
  diffvalues1 = gcn1[diffpairs]
  diffvalues2 = gcn2[diffpairs]
  diffgenes = cbind(diffgenes1,diffgenes2)
  diffgenes = t(apply(diffgenes, 1, sort))
  diff = cbind(diffgenes, diffvalues1, diffvalues2)
  diff = diff[!duplicated(diff[,c(1,2)]),c(1,2)]

  simipairs = which(diffnet<low & diffnet>0, arr.ind = TRUE)
  simipairs  = subset(simipairs, simipairs[,1]!=simipairs[,2])
  simigenes1 = genenames[simipairs[,1]]
  simigenes2 = genenames[simipairs[,2]]
  simivalues1 = gcn1[simipairs]
  simivalues2 = gcn2[simipairs]
  simigenes = cbind(simigenes1, simigenes2)
  simigenes = t(apply(simigenes, 1, sort))
  simi = cbind(simigenes, simivalues1, simivalues2)
  simi = simi[!duplicated(simi[,c(1,2)]),c(1,2)]

  return(list(similar = simi, different = diff))
}

#' @import ggplot2
#' @export
# plot the genes
scatter_genepairs <- function(dirname, extrainfo, cellcolors = NULL, datalist, genepairs_list, cellinfolist,
                              corr = c("bayescor", "pearson", "spearman", "cosine", "bigscale")){

  all = list()
  if(is.null(cellcolors))cellcolors =  brewer.pal(n = 9, name = "Set1")


  for(i in 1:nrow(genepairs_list)){
    genepair = genepairs_list[i,]
    for(j in 1:length(datalist)){
      data = datalist[[j]]
      data = log2(10^6*t(t(data)/colSums(data))+1)
      gc = gcn(data, norm = FALSE, CPM = FALSE, genepair, name = corr)
      data_genes = data[genepair,]
      data_corr = round(gc[1,2], 3)
      if(j==1){
        allgenes = data.frame(t(data_genes), newcelltype = cellinfolist[[j]]$newcelltype,
                              group = as.factor(c(rep(paste0(names(datalist)[j], " ", corr, "=", data_corr), ncol(data)))))
      }
      else allgenes = rbind(allgenes, data.frame(t(data_genes), newcelltype = cellinfolist[[j]]$newcelltype,
                                                    group = as.factor(c(rep(paste0(names(datalist)[j], " ", corr, "=", data_corr), ncol(data))))))

    }
    all[[i]] = allgenes
  }
  for(k in 1:length(all)){
    print(ggplot(all[[k]], aes(all[[k]][,1], all[[k]][,2], color = as.factor(newcelltype)))+
        geom_point(alpha=.4, size = .1) +
        xlab(colnames(all[[k]][1])) +
        ylab(colnames(all[[k]][2])) +
        theme(
          axis.text.y = element_text(colour="grey40", size=8, face="bold"),
          axis.text.x = element_text(colour="grey40", size=8, face="bold"),
          panel.background = element_rect(fill = "aliceblue",
                                        colour = "aliceblue"),
          panel.border = element_blank(),
          axis.title.x = element_text(colour="grey40", size=10, face="bold"),
          axis.title.y = element_text(colour="grey40", size=10, face="bold"))+
        scale_color_manual(values = cellcolors)+
        labs(color = 'Cell Type')+
        guides(color = guide_legend(override.aes = list(size = 3, alpha = .7)))+
      facet_grid(. ~ group))
    #ggsave(filename = paste0(dirname, "scatter", extrainfo, k, ".png"), width = 8, height = 3, units = "in")
  }

}

#' @import ggplot2
#' @export
histplot<-function(gcn1, gcn2, labels, xlab = "pearson correlation", main = "", yend = 7){
  dat = data.frame(dif = c(as.vector(gcn1), as.vector(gcn2)),
                   lab = as.factor(c(rep(1,nrow(gcn1)^2), rep(2,nrow(gcn2)^2))))
  print(ggplot(dat, aes(x=log2(dif+1), fill = lab)) +
          geom_histogram(alpha = 0.3, bins = 60, color = "gray30", position="identity", aes(y = ..density..)) +
          geom_segment(aes(x =0, y = 0, xend = 0, yend = yend), linetype="dashed",  color="black")+
          labs(x = xlab, y ="Density",title = main) +
          scale_fill_manual(labels = labels, values = c('orange', 'dodgerblue3'))+
          theme(legend.title=element_blank())+
          theme_light()+
          xlim(breaks = c(-1,1)))
}

#' #' @import ggplot2
#' #' @export
#' 
compare_simreal<-function(real, simlist, plotdir, realcellinfo = NULL, realdegeneinfo = NULL){
     
     # compare the group 
     if((!is.null(realcellinfo))&(!is.null(realdegeneinfo))){
       for(i in 1:(length(unique(colData(simlist[[1]])$Group)))){
         #  for marker genes in this cell type
         celltype = levels(realcellinfo$newcelltype)[i]
         subcells = as.character(realcellinfo$cell[which(realcellinfo$newcelltype==celltype)])
         
         subsimcells = lapply(simlist, function(sim)return(which(colData(sim)$Group==paste0("Group",i))))
         
         len = min(length(subcells), sapply(subsimcells, length))

         #  for all de genes in this cell type
         degenes = as.character(realdegeneinfo$genes)
         simdegenes = lapply(simlist, function(sim)return(rownames(assays(sim)$counts)[which(rowData(sim)$GeneGroup!="None")]))
         lende = min(length(degenes), sapply(simdegenes, length))
         
         subdereal = real[sample(degenes, lende), sample(subcells, len)]
         subdereal[1,which(colSums(subdereal)==0)] = 0.01
         
         subded =  list(Real = SingleCellExperiment(assays=list(counts=subdereal)))
         subdesimreal = lapply(as.list(1:length(simlist)), function(i){
           decounts = assays(simlist[[i]])$counts[sample(simdegenes[[i]], lende), sample(subsimcells[[i]], len)]
           decounts[1,which(colSums(decounts)==0)] = 0.01
           return(SingleCellExperiment(assays=list(counts=decounts)))
         })
         
         subded[2:(length(simlist)+1)] = subdesimreal
         names(subded)[2:(length(simlist)+1)] = names(simlist)

         compde <- compareSCEs(subded, point.size = 2, point.alpha = 0.2)
         diffde <- diffSCEs(subded, ref = "Real", point.size = 2, point.alpha = 0.2)
         resde <- list(Comp = compde, Diff = diffde)
         compareplot(resde, paste0(celltype,"_allde_zeisel"), plotdir)
         
         #  for all genes in this cell type
         
         suballreal = real[, sample(subcells, len)]
         suballreal[1,which(colSums(suballreal)==0)] = 0.01
         
         suballd =  list(Real = SingleCellExperiment(assays=list(counts=suballreal)))
         suballsimreal = lapply(as.list(1:length(simlist)), function(i){
                          decounts = assays(simlist[[i]])$counts[, sample(subsimcells[[i]], len)]
                          decounts[1,which(colSums(decounts)==0)] = 0.01
                          return(SingleCellExperiment(assays=list(counts=decounts)))
                        })
                        
         suballd[2:(length(simlist)+1)] = suballsimreal
         names(suballd)[2:(length(simlist)+1)] = names(simlist)

         compall <- compareSCEs(suballd, point.size = 2, point.alpha = 0.2)
         diffall <- diffSCEs(suballd, ref = "Real", point.size = 2, point.alpha = 0.2)
         resall <- list(Comp = compall, Diff = diffall)
         compareplot(resall, paste0(celltype,"_all_zeisel"), plotdir)
       }
      }
       
       suballallreal = real
       suballallreal[1,which(colSums(suballallreal)==0)] = 0.01
       suballalld =  list(Real = SingleCellExperiment(assays=list(counts=suballallreal)))
       suballallsimreal = lapply(as.list(1:length(simlist)), function(i){
                       decounts = assays(simlist[[i]])$counts
                       decounts[1,which(colSums(decounts)==0)] = 0.01
                       return(SingleCellExperiment(assays=list(counts=decounts)))
                       })
       suballalld[2:(length(simlist)+1)] = suballallsimreal
       names(suballalld)[2:(length(simlist)+1)] = names(simlist)
                       
       compallall <- compareSCEs(suballalld, point.size = 2, point.alpha = 0.2)
       diffallall <- diffSCEs(suballalld, ref = "Real", point.size = 2, point.alpha = 0.2)
       resallall <- list(Comp = compallall, Diff = diffallall)
       compareplot(resallall, paste0("zeisel"), plotdir)
     }
  



#'@import ggplot2 
#'@import cowplot
#'@export  
compareplot<-function(res, extrainfo, plotdir){
  plt1 <- makeCompPanel(res$Comp, title = extrainfo)
  plt1.name <- paste0("additional", extrainfo, "_comp")
  save_plot(paste0(plotdir, plt1.name, ".png"), plt1,
            ncol = 7, nrow = 10, base_height = 3)
  plt2 <- makeDiffPanel(res$Diff, title = extrainfo)
  plt2.name <- paste0("additional", extrainfo, "_diff")
  save_plot(paste0(plotdir, plt2.name, ".png"), plt2,
            ncol = 7, nrow = 10, base_height = 3)
  
  plots <- list(MeansComp   = res$Comp$Plots$Means,
                MeansDiff   = res$Diff$Plots$Means,
                VarsComp    = res$Comp$Plots$Variances,
                VarsDiff    = res$Diff$Plots$Variances,
                MeanVarComp = res$Comp$Plots$MeanVar,
                MeanVarDiff = res$Diff$Plots$MeanVar,
                LibSizeComp = res$Comp$Plots$LibrarySizes,
                LibSizeDiff = res$Diff$Plots$LibrarySizes)
  
  # Select the colours we are going to use
  cols <- c("black", "orange", "steelblue", "olivedrab")
  
  # Labels for datasets
  labels = c("Real" = "Real", "GroupSim" = "GroupSim", "HomoSim" = "HomoSim", "GeoSim = GeoSim")
  
  # Make adjustments to the plots
  for (idx in seq_along(plots)) {
    name <- names(plots)[idx]
    plot <- plots[[idx]]
    
    # Set a newnew theme
    plot <- plot +
      theme_cowplot(font_size = 12) +
      theme(legend.position = "none")
    
    # Set the colours, diff plots have one less dataset
    if (grepl("Comp", name)) {
      plot <- plot + scale_color_manual(values = cols, labels = labels) +
        scale_fill_manual(values = cols, labels = labels)
    } else {
      plot <- plot + scale_color_manual(values = cols[-1], labels = labels[-1]) +
        scale_fill_manual(values = cols[-1], labels = labels[-1])
    }
    
    # Boxplots are replotted with different properties, axis text adjusted and
    # x label removed
    if (!grepl("MeanVar", name)) {
      plot <- plot + geom_boxplot(aes(fill = Dataset),
                                  size = 1.5, alpha = 0.2) +
        scale_x_discrete(labels = labels[-1]) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    }
    
    plots[[idx]] <- plot
  }
  
  # Get a single legend to use
  leg <- get_legend(plots[["MeanVarComp"]] + theme(legend.position = "bottom"))
  
  # Assemble the panel
  panel <- ggdraw() +
    draw_label("A", 0.01, 0.986,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$MeansComp,  0.00, 0.77, 0.49, 0.23) +
    draw_label("B", 0.51, 0.986,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$MeansDiff,  0.51, 0.77, 0.49, 0.23) +
    draw_label("C", 0.01, 0.746,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$VarsComp,  0.00, 0.53, 0.49, 0.23) +
    draw_label("D", 0.51, 0.746,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$VarsDiff,  0.51, 0.53, 0.49, 0.23) +
    draw_label("E", 0.01, 0.506,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$MeanVarComp,  0.00, 0.29, 0.49, 0.23) +
    draw_label("F", 0.51, 0.506,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$MeanVarDiff,  0.51, 0.29, 0.49, 0.23) +
    draw_label("G", 0.01, 0.266,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$LibSizeComp,  0.00, 0.05, 0.49, 0.23) +
    draw_label("H", 0.51, 0.266,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$LibSizeDiff, 0.51, 0.05, 0.49, 0.23) +
    draw_plot(leg, 0.00, 0.00, 1.00, 0.04)
  
  save_plot(paste0(plotdir, extrainfo,  "_means.png"), panel, ncol = 2, nrow = 4)
  
  #====== zeros
  plots <- list(ZerosCellComp = res$Comp$Plots$ZerosCell,
                ZerosCellDiff = res$Diff$Plots$ZerosCell,
                ZerosGeneComp = res$Comp$Plots$ZerosGene,
                ZerosGeneDiff = res$Diff$Plots$ZerosGene,
                MeanZerosComp = res$Comp$Plots$MeanZeros,
                MeanZerosDiff = res$Diff$Plots$MeanZeros)
  
  # Make adjustments to the plots
  for (idx in seq_along(plots)) {
    name <- names(plots)[idx]
    plot <- plots[[idx]]
    
    # Set a newnew theme
    plot <- plot +
      theme_cowplot(font_size = 12) +
      theme(legend.position = "none")
    
    # Set the colours, diff plots have one less dataset
    if (grepl("Comp", name)) {
      plot <- plot + scale_color_manual(values = cols, labels = labels) +
        scale_fill_manual(values = cols, labels = labels)
    } else {
      plot <- plot + scale_color_manual(values = cols[-1], labels = labels[-1]) +
        scale_fill_manual(values = cols[-1], labels = labels[-1])
    }
    
    # Boxplots are replotted with different properties, axis text adjusted and
    # x label removed
    if (!grepl("MeanZeros", name)) {
      plot <- plot + geom_boxplot(aes(fill = Dataset),
                                  size = 1, alpha = 0.2) +
        scale_x_discrete(labels = labels[-1]) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    }
    
    plots[[idx]] <- plot
  }
  
  # Get a single legend to use
  leg <- get_legend(plots[["MeanZerosComp"]] + theme(legend.position = "bottom"))
  
  # Assemble the panel
  panel <- ggdraw() +
    draw_label("A", 0.01, 0.982,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$ZerosCellComp,  0.00, 0.69, 0.49, 0.31) +
    draw_label("B", 0.51, 0.982,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$ZerosCellDiff,  0.51, 0.69, 0.49, 0.31) +
    draw_label("C", 0.01, 0.662,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$ZerosGeneComp,  0.00, 0.37, 0.49, 0.31) +
    draw_label("D", 0.51, 0.662,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$ZerosGeneDiff,  0.51, 0.37, 0.49, 0.31) +
    draw_label("E", 0.01, 0.342,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$MeanZerosComp,  0.00, 0.05, 0.49, 0.31) +
    draw_label("F", 0.51, 0.342,
               fontface = "bold", hjust = 0, vjust = 0) +
    draw_plot(plots$MeanZerosDiff,  0.51, 0.05, 0.49, 0.31) +
    draw_plot(leg, 0.00, 0.00, 1.00, 0.04)
  
  save_plot(paste0(plotdir, extrainfo, "_zeros.png"), panel, ncol = 2, nrow = 3)
  
}

#'@import umap
#'@import ggplot2
#'@export 
umapplot<-function(data, celltype, labels, config = umap.defaults){
  embedding <- umap(t(data), config = config)
  dat = data.frame(UMAP1 = embedding$layout[,1], UMAP2 = embedding$layout[,2], celltype = celltype)
  cellcolors =  brewer.pal(n = 9, name = "Set1")
  
  print(ggplot(dat, aes(UMAP1, UMAP2, color = as.factor(celltype)))+
    geom_point(alpha=.6, size = .3) +
    theme(
      axis.text.y = element_text(colour="grey40", size=8, face="bold"),
      axis.text.x = element_text(colour="grey40", size=8, face="bold"),
      # panel.background = element_rect(fill = "aliceblue",
      #                                 colour = "aliceblue"),
      panel.border = element_blank(),
      axis.title.x = element_text(colour="grey40", size=10, face="bold"),
      axis.title.y = element_text(colour="grey40", size=10, face="bold"))+
    scale_color_manual(values = cellcolors, labels = labels)+ 
    labs(color = 'Cell Type')+
    guides(color = guide_legend(override.aes = list(size = 3, alpha = .7)))
  )
}  



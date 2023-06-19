## This R-file contains functions required for running the site-to-protein normalization script ("Normalize_MS2SiteToProtein.Rmd")




### PCA plot function ###
PCA_plot <- function(m, groups, samplenames, batch, legend_colors, plot_path=NULL, title=""){
  
  # load packages
  library(ggplot2)
  library(plotly)
  
  # replace NAs with 0
  m[is.na(m)] <- 0
  
  # calculate PCA
  pca_res <- prcomp(t(m))
  rot_mat <- pca_res$rotation
  res_final <- as.matrix(scale(t(m), center=TRUE, scale=FALSE)) %*% rot_mat
  eigenvalues <- pca_res$sdev^2
  anteil_var_pca1 <- round(eigenvalues[1]/sum(eigenvalues),digits=3)
  anteil_var_pca2 <- round(eigenvalues[2]/sum(eigenvalues),digits=3)
  
  ## create groups
  groups <- factor(groups, levels=names(legend_colors))
  colors <- legend_colors
  names(colors) <- levels(groups)
  colors <- colors[names(colors) %in% groups]
  
  ## ggplot PCR
  if(is.null(batch)){
    df_gg <- as.data.frame(res_final)
    df_gg$samplenames <- samplenames
    df_gg$groups <- groups
    gg <- ggplot(df_gg) + 
      geom_point(aes(x=PC1, y=PC2, col=groups, text=samplenames),size=5) +
      scale_color_manual(values=colors)+
      xlab(paste0("PC1  ","(",anteil_var_pca1*100,"%",")")) +
      ylab(paste0("PC2  ","(",anteil_var_pca2*100,"%",")")) +
      ggtitle(title) +
      theme_bw()
  } else {
    df_gg <- as.data.frame(res_final)
    df_gg$samplenames <- samplenames
    df_gg$groups <- groups
    df_gg$batch <- as.factor(batch)
    gg <- ggplot(df_gg) + 
      geom_point(aes(x=PC1, y=PC2, col=groups, text=samplenames, shape=batch),size=5) +
      scale_color_manual(values=colors)+
      xlab(paste0("PC1  ","(",anteil_var_pca1*100,"%",")")) +
      ylab(paste0("PC2  ","(",anteil_var_pca2*100,"%",")")) +
      ggtitle(title) +
      theme_bw()
  }
  
  # save plot
  if(!is.null(plot_path)){
    ggsave(plot=gg, filename=plot_path, width = 6, height = 4)
  }
  
  # print plot
  ggplotly(gg)
}



### Heatmap plot function
heatmap_plot <- function(m, groups, legend_colors, samplenames, type="normal", dendrogram="column", labrow="", bool_rowv=TRUE, bool_colv = TRUE, plot_path=NULL, title=""){
  
  # load packages
  library(gplots)
  
  # create groups
  names(colors) <- levels(groups)
  
  # replaces NAs with 0
  m[is.na(m)] <- 0
  colnames(m) <- samplenames
  
  # should rows be reordered
  if (bool_rowv){
    rowv <- as.dendrogram(hclust(dist(m)))
  } else {
    rowv <- FALSE
  }
  
  # reordering of columns
  if (is.logical(bool_colv)){
    if (bool_colv){
      colv <- as.dendrogram(hclust(dist(t(m))))
    } else {
      colv <- FALSE
    }
  } else {
    colv <- bool_colv
  }
  
  # specify colors
  if(is.null(legend_colors)){
    sidecolors <- rep("white", times=ncol(m))
  } else{
    sidecolors <- legend_colors[groups]
  }
  
  # create color palette
  heatmap_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  rdbu_colors = heatmap_pal(20)[2:19]
  heatmap_pal <- colorRampPalette(rdbu_colors)
  
  # plot heatmap if centered at 0
  if (type == "centered"){
    min_m <- min(m, na.rm=TRUE)
    max_m <- max(m, na.rm=TRUE)
    heatmap.2(m,         
              Rowv = rowv,
              Colv=colv,
              labRow=labrow, margins=c(8,8), ColSideColors = sidecolors, trace="none",col=heatmap_pal(50),
              breaks = seq(from=-2,to=2, length.out=51), 
              symkey = F,
              dendrogram=dendrogram, main=title)
    if(!is.null(plot_path)){
      pdf(file=plot_path, width = 6, height = 6)
      heatmap.2(m,         
                Rowv = rowv,
                Colv=colv,
                labRow=labrow, margins=c(8,8), ColSideColors = sidecolors, trace="none",col=heatmap_pal(50),
                breaks = seq(from=-2,to=2, length.out=51), 
                symkey = F,
                dendrogram=dendrogram,
                main=title)
      invisible(dev.off())
    }
  } 
  
  # plot heatmap if centered at 1
  if (type == "centered_at_1"){
    min_m <- min(m, na.rm=TRUE)
    max_m <- max(m, na.rm=TRUE)
    heatmap.2(m,         
              Rowv = rowv,
              Colv=colv,
              labRow=labrow, margins=c(8,8), ColSideColors = sidecolors, trace="none",col=heatmap_pal(50),
              breaks = seq(from=0,to=2, length.out=51), 
              symkey = F,
              dendrogram=dendrogram, main=title)
    if(!is.null(plot_path)){
      pdf(file=plot_path, width = 6, height = 6)
      heatmap.2(m,         
                Rowv = rowv,
                Colv=colv,
                labRow=labrow, margins=c(8,8), ColSideColors = sidecolors, trace="none",col=heatmap_pal(50),
                breaks = seq(from=0,to=2, length.out=51), 
                symkey = F,
                dendrogram=dendrogram,
                main=title)
      invisible(dev.off())
    }
  } 
  
  # plot heatmap if standardized
  if (type == "standardized"){
    min_m <- min(m, na.rm=TRUE)
    max_m <- max(m, na.rm=TRUE)
    heatmap.2(m,         
              Rowv = rowv,Colv= colv,labRow=labrow, margins=c(8,8), ColSideColors = sidecolors, trace="none",col=heatmap_pal(50),symkey = F,
              breaks = seq(from=-2,to=2, length.out=51),
              dendrogram=dendrogram, main=title)
    if (!is.null(plot_path)){
      pdf(file=plot_path, width = 6, height = 6)
      heatmap.2(m,         
                Rowv = rowv, Colv= colv,
                labRow=labrow, margins=c(8,8), ColSideColors = sidecolors, trace="none",col=heatmap_pal(50), symkey = F,
                breaks = seq(from=-2,to=2, length.out=51),
                dendrogram=dendrogram, main=title)
      invisible(dev.off())
    }
  }
}



### Write LOESS normalization function ###
loess_norm <- function(m){
  m[m==0] <- NA
  m_log <- log2(m)  
  m_norm <- 2^normalizeBetweenArrays(m_log, method="cyclicloess", cyclic.method = "fast")
  return(m_norm)
}



### Write DESeq normalization function (using DESeq2's size factor estimation) ###
DESeq_norm <- function(m, sizefactors=NULL){
  
  # create counts from intensity data in the required range
  m_copy <- m
  m_copy[is.na(m_copy)] <- 0
  m_counts <- round(log2(m_copy+1)*1000,digit=0)
  library(DESeq2)
  
  # if no sizefactors are supplied, calculate them based on m
  if (is.null(sizefactors)){
    # create an object summarized experiment class
    dds <- DESeqDataSetFromMatrix(countData = m_counts,
                                  colData = data.frame(condition=rep("group",times=ncol(m_counts))),
                                  design =  ~ 1)
    
    # calculate normalization factors via DESeq's estimateSizeFactors(). Save them in working directory (so they can be used later on a different table of the same experiment) 
    sizefactors <- estimateSizeFactors(dds)$sizeFactor
    if (!file.exists("Results")){
      dir.create("Results")
    }
    save(sizefactors, file=paste0(getwd(),"/Results/sizefactors.Rdata"))
  } else {
    sizefactors=sizefactors
  }
  
  # perform normalization by column-wise multiplication with size-factors
  m_counts_norm  <- sweep(m_counts, STATS=1/sizefactors, FUN="*", MARGIN = 2)
  
  # retransform to intensity range 
  m_norm <- 2^(m_counts_norm/1000) - 1
  
  # return normalized intensity matrix
  return(m_norm)
}




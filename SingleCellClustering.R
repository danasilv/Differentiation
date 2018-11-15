####Functions for the analysis of single-cell RNA-Seq data 2.#####
#UPDATED: August 22, 2018. KL.

#####################Load in R libraries#########################
library(Seurat)
library(Matrix)
library(stringr)
library(DiagrammeR)
#library(fifer)
library(NMF)
#library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamps)
#library(circlize)
library(ggplot2)
library(enrichR)
library(cluster)
library(tidyr)
#library(mygene)
library(dplyr)
library(igraph)
library(Rtsne)
library(scales)
#library(biomaRt)
library(diptest)
library(useful)
library(purrr)

#####################Load in dependencies#########################

#signatures = read.table("Resources/GBM_signatures.csv", header = TRUE, sep = ",", stringsAsFactors = F)
#signatures = as.list(signatures)
#signatures = lapply(signatures, function(x) x[!is.na(x)])
#gencode = read.table("Resources/gencode_v19_gene_pos.txt")


####################Single cell analysis via Seurat############################

#This function takes as inputs a tpm matrix and a boolean matrix of genetic information (derived from spiked-in primers). It converts the tpm to
#a seurat object, and runs the typical seurat preprocessing steps on that object. The genetic information is stored as metadata within the 
#seurat object.
#INPUTS: tpm contains the raw data. genetics contains the mutation information from amplified primers. hkgenes is the list of 
#housekeeping genes from Itay.

tpm.to.seurat = function(tpm, genetics, hkgenes = NULL, projectname = "Single Cell", figures_path = "results/",
                         mincells = 3, mingenes = 400)
{
        seurat_data <- CreateSeuratObject(raw.data = tpm, min.cells = mincells, min.genes = mingenes, project = projectname)
        genetics$well = gsub("-", "_", genetics$well)
        genetics = arrange(genetics, well)
        
        # Add in the genetic mutation data into the p6 seurat object, in metadata format. 
        genetics = genetics [genetics$well %in% intersect(genetics$well,seurat_data@cell.names),]
        
        # Create a dataframe that can be used to annotate the seurat object with the genetic information.
        genetics = full_join (data.frame(well = seurat_data@cell.names, stringsAsFactors = FALSE), genetics, all = TRUE)
        rownames(genetics) = genetics$well
        genetics = select(genetics, -well)
        
        # 0 = Not tested. 1 = absence of mutation. 2 = presence of mutation.
        genetics [!is.na(genetics)] = genetics [!is.na(genetics)] + 1
        genetics[is.na(genetics)] = 0
        
        seurat_data = AddMetaData(seurat_data, genetics)
        
        # Annotate the cells which contain genetic information within the seurat object.
        seurat_data = AddMetaData(seurat_data, data.frame(genetic.info = rowSums(genetics)!=0))
        
        # Filtering the seurat object, using mitochondrial genes and housekeeping genes.
        
        # Mitochondrial genes
        mito_genes <- grep(pattern = "^mt-", x = rownames(x = seurat_data@data), value = TRUE, ignore.case = T)
        percent_mito <- Matrix::colSums(seurat_data@raw.data[mito_genes, ]) / Matrix::colSums(seurat_data@raw.data)
        seurat_data <- AddMetaData(object = seurat_data, metadata = percent.mito, col.name = "percent.mito")
        
        # House keeping genes (list from Itay Tirosh)
        hkgenes <- as.vector(hkgenes$V1)
        hkgenes_found <- which(toupper(rownames(seurat_data@data)) %in% hkgenes)  # remove hkgenes that were not found
        n_expressed_hkgenes <- Matrix::colSums(seurat_data@data[hkgenes.found, ] > 0)
        seurat_data <- AddMetaData(object = seurat_data, metadata = n_expressed_hkgenes, col.name = "n.exp.hkgenes")
        
        # Output plots for filtering.
        VlnPlot(object = seurat_data, features.plot = c("nGene"))
        ggsave(paste0(figures_path, projectname, " nGene.png"), width = 10, height = 3)
        VlnPlot(object = seurat_data, features.plot = c("n.exp.hkgenes"))
        ggsave(paste0(figures_path, projectname, " nexphkgenes.png"), width = 10, height = 3)
        
        # Filtering.
        seurat_data <- FilterCells(object = seurat_data, subset.names = c("nGene", "n.exp.hkgenes"),
                                   low.thresholds = c(3000, 50), high.thresholds = c(Inf, Inf))
        
        # Normalization.
        seurat_data <- NormalizeData(object = seurat_data, normalization.method = "LogNormalize", scale.factor = 100000)
        
        # Find variable genes.
        seurat_data <- FindVariableGenes(seurat_data, mean.function = ExpMean, dispersion.function = LogVMR,
                                         x.low.cutoff = 0.1, x.high.cutoff = 7, y.cutoff = 1, do.plot = FALSE)
        
        p1 <- ggplot(seurat_data@hvg.info, aes(gene.mean, gene.dispersion)) + geom_point(size = 0.5, alpha = 1/10) + 
                labs(x = "average expression", y = "dispersion")
        p2 <- ggplot(seurat_data@hvg.info, aes(gene.mean, gene.dispersion.scaled)) + geom_point(size = 0.5, alpha = 1/10) + 
                labs(x = "average expression", y = "dispersion z-score")
        meanvar <- plot_grid(p1, p2, align = 'h', labels = c('A', 'B'))
        ggsave(paste0(figures_path, projectname, ' meanvar.png'), width = 5, height = 4)
        
        #Scaling data.
        seurat_data <- ScaleData(seurat_data, do.center = TRUE, do.scale = FALSE)
        
        #PCA
        seurat_data <- RunPCA(object = seurat_data, pc.genes = seurat_data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
        
        #Run Jackstraw
        seurat_data <- JackStraw(object = seurat_data, num.replicate = 100)
        
        #Output Jackstraw plot
        JackStrawPlot(object = seurat_data, PCs = 1:12)
        ggsave(paste0(figures_path, projectname, ' jackstraw.png'), width = 8, height = 8)
        
        return(seurat_data)
}

init <- function(seurat_data, hkgenes, projectname = "Single Cell", figures_path,
                    mincells = 3, mingenes = 400, minHkgenes = 50)
{
  
  # Mitochondrial genes
  mito_genes <- grep(pattern = "^mt-", x = rownames(x = seurat_data@data), value = TRUE, ignore.case = T)
  percent.mito <- Matrix::colSums(seurat_data@raw.data[mito_genes, ]) / Matrix::colSums(seurat_data@raw.data)
  seurat_data <- AddMetaData(object = seurat_data, metadata = percent.mito, col.name = "percent.mito")
  
  # House keeping genes (list from Itay Tirosh)
  hkgenes <- as.vector(hkgenes$V1)
  hkgenes.found <- which(toupper(rownames(seurat_data@data)) %in% hkgenes)  # remove hkgenes that were not found
  n_expressed_hkgenes <- Matrix::colSums(seurat_data@data[hkgenes.found, ] > 0)
  seurat_data <- AddMetaData(object = seurat_data, metadata = n_expressed_hkgenes, col.name = "n.exp.hkgenes")
  
  # Output plots for filtering.
  VlnPlot(object = seurat_data, features.plot = c("nGene"))
  ggsave(paste0(figures_path, projectname, "_nGene.png"), width = 10, height = 3)
  VlnPlot(object = seurat_data, features.plot = c("n.exp.hkgenes"))
  ggsave(paste0(figures_path, projectname, "_nexphkgenes.png"), width = 10, height = 3)
  
  # Filtering.
  seurat_data <- FilterCells(object = seurat_data, subset.names = c("nGene", "n.exp.hkgenes"),
                             low.thresholds = c(mingenes, minHkgenes), high.thresholds = c(Inf, Inf))
  
  # Normalization.
  seurat_data <- NormalizeData(object = seurat_data, normalization.method = "LogNormalize", scale.factor = 100000)
  
  
  return (seurat_data)
}

cluster <- function(seurat_data, hkgenes, projectname = "Single Cell", figures_path) {
  
  # Find variable genes.
  seurat_data <- FindVariableGenes(seurat_data, mean.function = ExpMean, dispersion.function = LogVMR,
                                   x.low.cutoff = 0.1, x.high.cutoff = 7, y.cutoff = 1, do.plot = FALSE)
  
  
  #Scaling data.
  seurat_data <- ScaleData(seurat_data, do.center = TRUE, do.scale = FALSE)
  
  #PCA
  seurat_data <- RunPCA(object = seurat_data, pc.genes = seurat_data@var.genes, do.print = TRUE, pcs.print = 1:15, genes.print = 5)
  
  #Run Jackstraw
  seurat_data <- JackStraw(object = seurat_data, num.replicate = 100)
  
  #Output Jackstraw plot
  JackStrawPlot(object = seurat_data, PCs = 1:15)
  ggsave(paste0(figures_path, projectname, '_jackstraw.png'), width = 8, height = 8)
  
  PCElbowPlot(object = seurat_data)
  ggsave(paste0(figures_path, projectname, '_elbow.png'), width = 8, height = 8)
  
  PCHeatmap(object = seurat_data, pc.use = 1:15, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  ggsave(paste0(figures_path, projectname, '_PCHeatMap.png'), width = 8, height = 8)
  
  return(seurat_data)
}       

#This function takes as an argument a seurat object that has been filtered, and returns calculated scores for a list of genes For control genes, we are binning all genes into 25 bins of expression levels, and randomly
#selecting 100 genes from the same expression bin as each signature gene (see Tirosh et. al., Nature). The output of the function is a dataframe
#of signature scores that may be used to annotate a seurat object as metadata.

calc.sig.genes = function(seurat_data, signature.genes, name)
{ 
  #Binning all detected genes into 25 bins to identify control genes for each signature.
  genes_detected <- rownames(seurat_data@data)
  aggregate_exprs = rowMeans(as.matrix(seurat_data@data))
  
  # Kevin version: (somethomes it's mpossib;e to break the quantiles)
  #control_bins = cut(aggregate_exprs, breaks = c(quantile(aggregate_exprs, probs = seq(0,1,by=0.04))),
  #                   labels = 1:25, include.lowest=TRUE)
  # Dana changed:
  control_bins = cut(aggregate_exprs, breaks = c(seq(floor(min(aggregate_exprs)), ceiling(max(aggregate_exprs)), length.out = 25)),
                     labels = 1:24, include.lowest=TRUE)
  control_bins = data.frame(genes = genes_detected, mean = aggregate_exprs, bin = control_bins, stringsAsFactors = FALSE)
  
  binned_genes = list()
  for (i in 1:25)
  {
    binned_genes[[i]] <- rownames(control_bins)[control_bins$bin == i]
  }
  
  newsigs_genes_exprs <- FetchData(seurat_data, signature.genes)
  
  newsigs_mean <- rowMeans(newsigs_genes_exprs)
  
  
  #For the gene signature, obtain the control genes and calculate the gene signature for each cell.
 
  # Dana: ToDo - merge this with gbm2.
  temp = control_bins[control_bins$genes %in% signature.genes,]
  #For each gene on the list, pick out 100 genes randomly from that bin.
  temp_control_genes = vector()
  for (j in 1:nrow(temp))
  {
    temp2 = sample(binned_genes[[temp[j,]$bin]],min(100,length(binned_genes[[temp[j,]$bin]])))
    temp_control_genes = c(temp_control_genes, temp2)
  }
  # For each cell, want the mean expression level across all the control genes.
  temp_control_genes = FetchData(seurat_data, temp_control_genes)
  control_genes_exprs = rowMeans(temp_control_genes)
  # Subtract control gene expression from signature gene expression to get scores.
  newsigs_mean = newsigs_mean - control_genes_exprs
 
  
  #format the dataframe for output
  newsigs_mean = data.frame(newsigs_mean)
  colnames(newsigs_mean) <- name
  
  return(newsigs_mean)
  
}

#This function takes as an argument a seurat object that has been filtered, and returns calculated scores for the new signatures
#from GBM2.0, and the old Verhaak signatures. For control genes, we are binning all genes into 25 bins of expression levels, and randomly
#selecting 100 genes from the same expression bin as each signature gene (see Tirosh et. al., Nature). The output of the function is a dataframe
#of signature scores that may be used to annotate a seurat object as metadata.

gbm2.sigs = function(seurat_data, signature.dir ="/Volumes/ahg_regevdata2/projects/Glioma_scGenetics/resources/genesignatures/")
{ 
        #Binning all detected genes into 25 bins to identify control genes for each signature.
        genes_detected <- rownames(seurat_data@data)
        aggregate_exprs = rowMeans(as.matrix(seurat_data@data))
        
        # Taking the mean expression of each gene across the data set, and binning the genes based on their mean expression
        control_bins = cut(aggregate_exprs, breaks = c(quantile(aggregate_exprs, probs = seq(0,1,by=0.04), na.rm = TRUE)),
                                             labels = 1:25, include.lowest=TRUE, na.rm = TRUE)
          #cut(aggregate_exprs, breaks = c(seq(floor(min(aggregate_exprs)), ceiling(max(aggregate_exprs)), length.out = 25)),
          #                 labels = 1:24, include.lowest=TRUE)
        control_bins = data.frame(genes = genes_detected, mean = aggregate_exprs, bin = control_bins, stringsAsFactors = FALSE)
        
        binned_genes = list()
        for (i in 1:25)
        {
                binned_genes[[i]] <- rownames(control_bins)[control_bins$bin == i]
        }

        #A. Look at new stemness markers as defined by GBM2.0
        ACvMes <- read.csv(paste0(signature.dir,"ACvMES.csv"), skip=1)
        colnames(ACvMes)[3] <- "mes_ac"
        NPCvOPC <- read.csv(paste0(signature.dir,"NPCvOPC.csv"), skip=1)
        colnames(NPCvOPC)[3] <- "npc_opc"
        
        #reformatting. 
        newsigs <- list()
        for (i in 1:3) {newsigs[[i]] <- ACvMes[,i]}
        for (i in 1:3) {newsigs[[i+3]] <- NPCvOPC[,i]}
        for (i in 1:6) {newsigs[[i]] <- newsigs[[i]][newsigs[[i]] %in% genes_detected]}
        
        newsigs_genes_exprs <- list()
        for (i in 1:6) {newsigs_genes_exprs[[i]] <- FetchData(seurat_data, newsigs[[i]])}
        
        # For ech gene in the signature -= copmpute its mean expression
        newsigs_mean <- list()
        for (i in 1:6) {newsigs_mean[[i]] <- rowMeans(newsigs_genes_exprs[[i]])}
        
        
        #For each gene signature, obtain the control genes and calculate the gene signature for each cell.
        for (i in 1:length(newsigs))
        {
               
                # Todo: check what happens when we have a gene in the signature that was not detected
                temp = control_bins[control_bins$genes %in% newsigs[[i]],]
                
                #For each gene on the list, pick out 100 genes randomly from that bin.
                temp_control_genes = vector()
                for (j in 1:nrow(temp))
                {
                        temp2 = sample(binned_genes[[temp[j,]$bin]],100)
                        temp_control_genes = c(temp_control_genes, temp2)
                }
                
                # For each cell, want the mean expression level across all the control genes.
                temp_control_genes = FetchData(seurat_data, temp_control_genes)
                control_genes_exprs = rowMeans(temp_control_genes)

                # Subtract control gene expression from signature gene expression to get scores.
                newsigs_mean[[i]] = newsigs_mean[[i]] - control_genes_exprs
        }
        
        #format the dataframe for output
        newsigs_mean = data.frame(newsigs_mean)
        colnames(newsigs_mean) =  c(colnames(ACvMes), colnames(NPCvOPC))
        
        return(newsigs_mean)
}

####################Single cell analysis via hierarchical clustering/NMF##############################
#
#This function converts a .merged.txt file to the tpm format that will be used in subsequent steps.
convert.merged.file = function(data)
{
        data = dplyr::select(data, -X2transcript_id)
        rownames(data) = data$X1gene_id
        data = dplyr::select(data, -X1gene_id)
        return(data)
}


#This function takes as argument a tpm matrix and a vector of Sample IDs (sample_ident).
#It filters the data, and centers its expression levels.
#The function also generates two QC plots. One for the number of genes expressed, one for aggregate expression of genes.
#The output is a list containing filtered tpm matrix ($TPM) and the sample identities of the filtered tpm matrix ($sample_ident).
tpm.process = function(tpm, sample_ident, plot_path = "figures/", nGene_cutoff_low = 2500, nGene_cutoff_high = 7000, Ea_cutoff = 4)
{
        #the "temp" is for the second step of filtering (aggregate gene expression levels)
        temp = tpm
        
        #Convert expression levels.
        tpm = log2(tpm/10 + 1)
        
        #1. Filter cells based on the number of genes expressed.
        nGene = apply(tpm, 2, function(x) sum(x != 0))
        
        #Generate QC plots
        nGene_toplot = data.frame(Sample = sample_ident, nGene = nGene)
        nGene_plot = ggplot(nGene_toplot, aes(factor(Sample), nGene)) + geom_violin() + geom_jitter()
        nGene_plot = nGene_plot + xlab("Sample") + ylab("Number of genes expressed")
        nGene_plot = nGene_plot + theme(axis.text.x = element_text(angle = 60, hjust = 1))
        nGene_plot = nGene_plot + geom_hline(yintercept = nGene_cutoff_low, linetype = "dashed", color = "red")
        nGene_plot = nGene_plot + geom_hline(yintercept = nGene_cutoff_high, linetype = "dashed", color = "red")
        print(nGene_plot)
        ggsave(paste0(plot_path, "nGene QC.pdf"), height = 5.4, width = 9)
        
        keep = (nGene >= nGene_cutoff_low) & (nGene <= nGene_cutoff_high)
        tpm = tpm[,keep]
        sample_ident = sample_ident[keep]
        
        #2. Filter genes based on aggregate expression.
        Ea = apply(temp, 1, function(x) log2(mean(x)+1))
        
        #   Generate QC plot.
        Ea_toplot = data.frame(Aggregate_Expression = Ea)
        Ea_plot = ggplot(Ea_toplot, aes(Aggregate_Expression)) + geom_histogram(binwidth = 0.1)
        Ea_plot = Ea_plot + xlab("Aggregate Gene Expression") + ylab("Count")
        Ea_plot = Ea_plot + geom_vline(xintercept = Ea_cutoff, linetype = "dashed", color = "red")
        print(Ea_plot)
        ggsave(paste0(plot_path, "Ea QC.pdf"), height = 5.4, width = 9)
        
        keep_gene = Ea >= Ea_cutoff
        tpm = tpm[keep_gene,]
        rm(temp)
        
        #Center gene expression.
        mean_exprs = apply(tpm, 1, mean) 
        tpm = sweep(tpm, 1, mean_exprs, "-")
        
        #Return results
        tpm_return = list()
        tpm_return[[1]] = tpm
        tpm_return[[2]] = sample_ident
        names(tpm_return) = c("TPM", "sample_ident")
        return (tpm_return)
}

#Code is adapted from Chris Rodman, Summer 2018.
#This function will take as argument a filtered tpm matrix (from tpm.process), perform hierarchical clustering, and 
#generate plots (heatmap of pearson correlations, and hierachical clustering) for exploratory analysis of the data.
#if marker genes are supplied, a separate plot will be generated at the bottom of the heatmap.
#It will return the hclust object.
# Input: 
# tpm - logged scale TPM
# sample_ident - vector of sample names, for example MGH143
# colours - a vector of color names, its length is the number of unique samples in the tpm
# marker_genes - if specified will annotate the cells with marker genes in the bottom
tpm.cluster = function(tpm, sample_ident, plot_path = "figures/", colours,
                       marker_genes = 0)
{ 
        C = cor(tpm)
        hc = hclust(as.dist(1 - C), method = "average")
        
        #Annotate correlation heatmap by Sample.
        ha = rowAnnotation(df = data.frame(Sample = sample_ident), 
                               show_annotation_name = TRUE,
                               col = colours,
                               annotation_name_gp = gpar(fontsize = 8),
                               annotation_height = unit.c(unit(0.75, "cm"), unit(0.75, "cm")))
        
        #Remove marker genes that are not found in the data.
        if (marker_genes[1] != 0)
        {
                marker_genes = intersect(marker_genes, rownames(tpm))
        }
        
        #Generate heatmap
        heatmap = Heatmap(C, 
                col = colorRamp2(seq(-0.05, 0.15, length.out=299), matlab.like2(299)),
                name = "Pearson", 
                show_row_dend = FALSE, 
                show_column_dend = FALSE,
                row_dend_reorder = FALSE, 
                column_dend_reorder = FALSE, 
                cluster_rows = hc, 
                cluster_columns = hc, 
                show_row_names = FALSE, 
                show_column_names = FALSE, 
                width = 6,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"), 
                                            labels_gp = gpar(fontsize = 8)))
        
        #Generation heatmap of marker genes.
        if (marker_genes[1] != 0)
        { 
                marker_genes_exprs = tpm[marker_genes,]
                heatmap2 = Heatmap(t(marker_genes_exprs),
                          col = colorRamp2(seq(-6, 8, length.out=299), rev(colorRampPalette(brewer.pal(11, "RdBu"))(299))),
                          cluster_rows = FALSE,
                          show_row_dend = FALSE, 
                          show_column_dend = FALSE,
                          row_dend_reorder = FALSE, 
                          column_dend_reorder = FALSE,
                          cluster_columns = FALSE,
                          show_column_names = TRUE,
                          show_row_names = FALSE,
                          column_names_gp = gpar(fontsize = 12),
                          width = 1
                          )
                png(paste0(plot_path, "Correlation Heatmap.png"), width = 900, height = 700)
                print(heatmap + heatmap2 + ha)
                dev.off()
        }else{ 
        png(paste0(plot_path, "Correlation Heatmap.png"), width = 900, height = 700)
        print(heatmap + ha)
        dev.off()
        }
        
        #Plot Clustering on its own
        pdf(paste0(plot_path, "Hierarchical Clustering.pdf"), width = 14, height = 7)
        hierarchical_clustering = plot(as.hclust(cut(as.dendrogram(hc),0.85)$upper), labels = FALSE, 
                                                  ylab = "Pearson Correlation",
                                       xlab = FALSE)
        print(hierarchical_clustering)
        dev.off()
        
        #Plot heatmap of marker genes
        if (marker_genes[1] != 0)
        { 
                pdf(paste0(plot_path, "Marker Genes.pdf"), width = 7, height = 14)
                print(hierarchical_clustering + heatmap2)
                dev.off()
        }
        
        return(hc)
}

#Code is adapted from Chris Rodman, Summer 2018.
#This function will take as argument a filtered tpm matrix (from tpm.process), a clustering object (from tpm.cluster),
#the number of desired clusters, and return a dataframe containing the genes that are most highly expressed by cells within
#each of the clusters.
# Input: K- number of clusters
# n_gene number of genes to return
tpm.expressed.genes = function(tpm, hc, k = 3, n_gene = 50)
{
        hc_clusters = cutree(hc, k = k)
        
        top_genes = list()
        for (i in 1:k)
        { 
                top_genes[[i]] = names(sort(apply(tpm[ ,hc_clusters == i], 1, mean), decreasing = TRUE)[1:n_gene])
        }
        top_genes = data.frame(top_genes)
        colnames(top_genes) = paste0("Cluster ", seq(1:k))
        
        return(top_genes)
}

vioplot2<-function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                    horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                    lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                    at, add = FALSE, wex = 1, drawRect = TRUE) 
{
  if(!is.list(x)){
    datas <- list(x, ...)
  } else{
    datas<-x
  }
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

calc.nmf = function(df.data,name,results.dir,num.nmf.factors = 4,num.nmf.top.genes = 30) {
  t.df.TimePoint <- t(df.data)
  t.df.TimePoint <- t.df.TimePoint[rowSums(t.df.TimePoint) > 0,]
  nmf.df.TimePoint <- nmf(as.matrix(t.df.TimePoint),num.nmf.factors)
  w = basis(nmf.df.TimePoint)
  signatures = apply(w, 2, function(x) head(rownames(w)[order(x, decreasing = TRUE)],num.nmf.top.genes))
  
  write.table(signatures,paste0(results.dir,"/nmf_",name,"_",num.nmf.factors,"factors_top_",num.nmf.top.genes,"genes.txt"), row.names = FALSE)
}

prepare.data.for.breaks = function(vec.expression, dims)
{
  return(log(vec.expression + 1) + rnorm(dims, 0, 1e-6))
}


# Creates a violin plot for each gene in the lists, for each data frame in df.TimePoints
# input:
# df.TimePoints - list of data frames. Each data frame is a tpm matrix
# names - list of names corresponding to the data frames
plot.volin.for.TF = function(df.TimePoints,names.TimePoints, filename, results.dir,Suva_genes = c("OLIG1","OLIG2","SOX2","SOX8","POU3F3","HES6","POU3F2","HEY2","SOX5","RFX4","KLF15","CITED1","LHX2","VAX2","MYCL","SALL2","HES1","SOX4","DLX2","ID4","ID2","MYC","HEY1","HEY1")
                             ,Filbin_stemness = c("CD24","SOX6","SOX10")
                             ,Filbin_differentiation = c("IFI6","IFI16","IFITM2","STAT2","S100A4","S100A6","S100A10","S100A11","S100A13","S100A16","SERPING1","SERPINGH1","SERPINGI1"))
{
  #Suva_genes = c("OLIG1","OLIG2","SOX2","SOX8","POU3F3","HES6","POU3F2","HEY2","SOX5","RFX4","KLF15","CITED1","LHX2","VAX2","MYCL","SALL2","HES1","SOX4","DLX2","ID4","ID2","MYC","HEY1","HEY1")
  #Filbin_stemness = c("CD24","SOX6","SOX10")
  #Filbin_differentiation = c("IFI6","IFI16","IFITM2","STAT2","S100A4","S100A6","S100A10","S100A11","S100A13","S100A16","SERPING1","SERPINGH1","SERPINGI1")
  Filbin_genes = c(Filbin_stemness,Filbin_differentiation)
  
  genes<-Filbin_genes  

  
# for (df.TimePoint in df.TimePoints) {
#    for (gene in c(union(Filbin_genes, Suva_genes))) {
#      if (gene %in% colnames(df.TimePoints[[1]])) {
#        print (gene)
#        df.TimePoint[,gene] <- log(df.TimePoint[,gene] + 1) + rnorm(dim(df.TimePoint)[1], 0, 1e-6)
#      }
#    }
#  }
  
  
  pdf (paste0(results.dir,filename,"_Filbin.pdf"),width = 20, height = 20)
  layout(matrix(c(1:length(genes)), nr=divisors(length(genes))[3], byrow=F))
  
  # For each gene create the violin plots, DIPG needs to be logged scaled
  #ToDo: calc the ylim auto
  for (gene in genes) 
  {
    print(gene)
    if (gene %in% colnames(df.TimePoints[[1]])) {

       if (gene %in% Filbin_stemness) {
        
        vioplot2(lapply(df.TimePoints, function(x) prepare.data.for.breaks(x[,gene], length(x[,gene]))), names=names.TimePoints, col = c("gold"), ylim = c(0,1.5))
      }
      else {
        vioplot2(lapply(df.TimePoints, function(x) prepare.data.for.breaks(x[,gene], length(x[,gene]))), names=names.TimePoints, col = c("green"), ylim = c(0,1.5))
      }
      title(ylab=gene)
    }
  }
  dev.off()
  
  pdf (paste0(results.dir,filename,"_Suva.pdf"),width = 20, height = 20)
  layout(matrix(c(1:length(genes)), nr=4, byrow=F))
  
  # For each gene create the violin plots, DIPG needs to be logged scaled
  for (gene in Suva_genes) 
  {
    print(gene)
    if (gene %in% colnames(df.TimePoints[[1]])) {
      
      vioplot2(lapply(df.TimePoints, function(x) prepare.data.for.breaks(x[,gene], length(x[,gene]))), names=names.TimePoints, col = c("light blue"))
      title(ylab=gene)
    }
  }
  dev.off()
}

#This function takes as argument a filtered tpm matrix (from tpm.process, or after removing non-malignant cells), and then 
#runs nmf to identify gene signatures. The output is a list containing gene names for each obtained signature.
#Input: filtered tpm matrix. Sample identity vector.
# n_gene_per_signature: the number of genes to return per signature
tpm.to.nmf = function(tpm_filtered, sample_ident, n_gene_per_signature = 30)
{
        #convert negative values to 0, because NMF doesn't take negatives also we are more interested in positively expressed genes.
        tpm_filtered[tpm_filtered<0] <- 0
        
        #Get vector of unique sample IDs
        sample = sample_ident
        unique_samples = unique(sample)
        
        #convert data to list form.
        tpm_filtered_listed = list()
        for (i in 1:length(unique_samples))
        {
                tpm_filtered_listed[[i]] = tpm_filtered[,sample == unique_samples[i]]
                tpm_filtered_listed[[i]] = tpm_filtered_listed[[i]][rowSums(tpm_filtered_listed[[i]]) > 0,]
        }
        
        #obtain gene signatures via NMF.
        signatures = list()
        for (i in 1:length(unique_samples))
        {
                print(paste0("Running NMF on Sample ", unique_samples[i]))
                temp = nmf (tpm_filtered_listed[[i]],10)
                w = basis(temp)
                signatures [[i]] = apply(w, 2, function(x) head(rownames(w)[order(x, decreasing = TRUE)],n_gene_per_signature))
        }
        
        names(signatures) = unique_samples
        return (signatures)        
}

#This function takes as argument a filtered tpm matrix (from tpm.process, or after removing non-malignant cells), and NMF
#signatures (from tpm.to.nmf), and calculates signature scores for each cell in the tpm matrix and outputs it. The method 
#for obtaining control genes is the same as the new.sigs function (see Tirosh et. al., Nature). 
calc.signature.scores = function(tpm_filtered, nmf_signatures)
{ 
        #Aggregate all signatures together, if present in list format.
        signatures_aggregate = data.frame(nmf_signatures)
        colnames(signatures_aggregate) = paste0("Signature",1:ncol(signatures_aggregate))
        
        #Calculate aggregate expression levels of all genes and bin by expression levels.
        genes = rownames(tpm_filtered)
        aggregate_exprs = rowMeans(tpm_filtered)
        control_bins = cut(aggregate_exprs, breaks = c(quantile(aggregate_exprs, probs = seq(0,1,by=0.04))),
                           labels = 1:25, include.lowest=TRUE)
        control_bins = data.frame(genes = genes, mean = aggregate_exprs, bin = control_bins, stringsAsFactors = FALSE)
        
        binned_genes = list()
        for (i in 1:25)
        {
                binned_genes[[i]] <- rownames(control_bins)[control_bins$bin == 1]
        }
        
        #Calculate a signature score for each cell.
        tpm_filtered = data.frame(tpm_filtered, gene = rownames(tpm_filtered), stringsAsFactors = FALSE)
        
        signature_scores  = list()
        
        for (i in 1:ncol(signatures_aggregate))
        {
                sig = signatures_aggregate[,i]
                test = inner_join(data.frame (gene = sig, stringsAsFactors = FALSE), tpm_filtered)
                test = colMeans(test[,2:ncol(test)])
                
                #Find the bins of each gene within each gene signature.
                temp <- inner_join(data.frame(genes = as.character(signatures_aggregate[,i]), stringsAsFactors = FALSE), control_bins)
                
                #For each gene on the list, pick out 100 genes randomly from that bin.
                temp_control_genes = vector()
                for (j in 1:nrow(temp))
                {
                        temp2 = sample(binned_genes[[temp[j,]$bin]],100)
                        temp_control_genes = c(temp_control_genes, temp2)
                }
                
                # For each cell, want the mean expression level across all the control genes.
                temp_control_genes = data.frame(index = 1:(100*length(sig)), gene = temp_control_genes)
                temp_control_genes = inner_join(temp_control_genes, tpm_filtered, by = "gene")
                
                control_genes_exprs = colMeans(temp_control_genes[,3:ncol(temp_control_genes)])
                
                #Calculate gene signature for each cell.
                test = test - control_genes_exprs
                signature_scores[[i]] <- test
        }
        tpm_filtered = dplyr::select(tpm_filtered, -gene)
        
        # Make a dataframe of all signature scores.
        signature_scores = data.frame(signature_scores)
        colnames(signature_scores) = paste0("Signature",1:ncol(signatures_aggregate))
        rownames(signature_scores) = colnames(tpm_filtered)
        
        return (signature_scores)
}

#This function is the same as the calc.signature.scores function, except it works with gene signatures that are supplied in a list
#format (i.e each signature can have a different number of genes). This is for example useful for the GBM 2.0 signatures, which
#are of varying length.
calc.signature.scores.list = function(tpm_filtered, signatures)
{ 
        
        #Calculate aggregate expression levels of all genes and bin by expression levels.
        genes = rownames(tpm_filtered)
        aggregate_exprs = rowMeans(tpm_filtered)
        control_bins = cut(aggregate_exprs, breaks = c(quantile(aggregate_exprs, probs = seq(0,1,by=0.04))),
                           labels = 1:25, include.lowest=TRUE)
        control_bins = data.frame(genes = genes, mean = aggregate_exprs, bin = control_bins, stringsAsFactors = FALSE)
        
        binned_genes = list()
        for (i in 1:25)
        {
                binned_genes[[i]] <- rownames(control_bins)[control_bins$bin == 1]
        }
        
        #Calculate a signature score for each cell.
        tpm_filtered = data.frame(tpm_filtered, gene = rownames(tpm_filtered), stringsAsFactors = FALSE)
        
        signature_scores  = list()
        
        for (i in 1:length(signatures))
        {
                sig = signatures[[i]]
                test = inner_join(data.frame (gene = sig, stringsAsFactors = FALSE), tpm_filtered)
                test = colMeans(test[,2:ncol(test)])
                
                #Find the bins of each gene within each gene signature.
                temp <- inner_join(data.frame(genes = as.character(sig, stringsAsFactors = FALSE)), control_bins)
                
                #For each gene on the list, pick out 100 genes randomly from that bin.
                temp_control_genes = vector()
                for (j in 1:nrow(temp))
                {
                        temp2 = sample(binned_genes[[temp[j,]$bin]],100)
                        temp_control_genes = c(temp_control_genes, temp2)
                }
                
                # For each cell, want the mean expression level across all the control genes.
                temp_control_genes = data.frame(index = 1:length(temp_control_genes), gene = temp_control_genes)
                temp_control_genes = inner_join(temp_control_genes, tpm_filtered, by = "gene")
                
                control_genes_exprs = colMeans(temp_control_genes[,3:ncol(temp_control_genes)])
                
                #Calculate gene signature for each cell.
                test = test - control_genes_exprs
                signature_scores[[i]] <- test
        }
        tpm_filtered = dplyr::select(tpm_filtered, -gene)
        
        # Make a dataframe of all signature scores.
        signature_scores = data.frame(signature_scores)
        colnames(signature_scores) = names(signatures)
        rownames(signature_scores) = colnames(tpm_filtered)
        
        return (signature_scores)
}

#This function simply visualizes signature scores (from calc.signature.scores) as a heatmap, and as a hierarchical clustering
#of signatures, for the purposes of identifying signatures to merge as meta-signatures. It outputs a clustered object of the signatures
#and generates plots of the heatmap and hierarchical clustering. The purpose of this function is to assist in manually combining
#signatures into meta-signatures using the create.metasignature function.
cluster.signature.scores = function(nmf_signature_scores, plot_path)
{
        signature_cluster = hclust(dist(t(nmf_signature_scores)))
        #Visualize.
        pdf (paste0(plot_path, "NMF signatures heatmap.pdf"), width = 8, height = 8)
        plot = Heatmap(as.matrix(nmf_signature_scores),
                cluster_rows = T,
                cluster_columns = T,
                show_row_names = F,
                column_names_gp = gpar(fontsize = 8),
                col = colorRamp2(seq(-3, 8, length.out=299), rev(colorRampPalette(brewer.pal(11, "RdBu"))(299)))
                )
        print(plot)
        dev.off()
        
        pdf (paste0(plot_path, "NMF Signature clusters.pdf"), width = 8, height = 4)
        plot(signature_cluster)
        dev.off()
        
        return(signature_cluster)
}

#This function takes as arguments the filtered tpm (from tpm.process), signatures (from tpm.to.nmf), signature scores (from calc.signature.scores), a boolean 
#index vector specifying which signatures comprise the desired meta-signature (cluster_indices), the desired name of the meta signature, and the desired number
#of genes for the meta signature. The cluster_indices is specifified manually, based on visual inspect of the heatmap (cluster.signature.scores).
#The function determines the genes to keep within the meta-signature by ordering of correlation with average cell score of the
#meta-signature. It will write out the genes of that meta signature, the signatures that comprise that meta signature, and create a plot for that meta signature.
#The output of the function is a data frame, containing the expression level of each cell for each gene of the meta-signature.

create.metasignature = function(tpm_filtered, plot_path, results_path, colours, nmf_signatures, signature_scores, cluster_indices, cluster_name, n_metasignature_genes = 30)
{ 
        #Aggregate all signatures together.
        signatures_aggregate = data.frame(nmf_signatures)
        colnames(signatures_aggregate) = paste0("Signature",1:ncol(signatures_aggregate))
        #write out a list of signature that comprise the meta signature.
        write.csv(colnames(signatures_aggregate)[cluster_indices], paste0(results_path, "MetaSignature ", cluster_name, " members.csv"))

        meta_signature = as.character(unlist(as.list(signatures_aggregate[,cluster_indices])))
        meta_signature = unique(meta_signature)
        #average cell scores of programs in the set.
        avg_scores = rowMeans(signature_scores[,cluster_indices])

        #ranked genes by their correlation with the average cell scores of the programs
        ranked = sapply (meta_signature, function(x) cor(as.numeric(tpm_filtered[which(rownames(tpm_filtered) == x),]),avg_scores))

        #Sort and keep the top genes for the meta signature
        keep_meta_signature = sort(ranked, decreasing = TRUE) [1:n_metasignature_genes]
        #Write out the genes belonging to that meta signature.
        write.csv(names(keep_meta_signature), paste0(results_path, "MetaSignature ", cluster_name, " genes.csv"))
        
        to_plot = as.matrix(tpm_filtered[rownames(tpm_filtered) %in% names(keep_meta_signature),]) 
        to_plot = to_plot[,order(colSums(to_plot))]
        
        #Get sample identities, for plotting purposes
        sample_ident = colnames(to_plot)
        sample_ident = unlist(lapply(strsplit(sample_ident, "_"), function(x) x[1]))
        
        keep_meta_signature = names(keep_meta_signature)
        
        #Annotate correlation heatmap by Sample.
        ha = rowAnnotation(df = data.frame(Sample = sample_ident), 
                           show_annotation_name = TRUE,
                           col = colours,
                           annotation_name_gp = gpar(fontsize = 6),
                           annotation_height = unit.c(unit(0.75, "cm"), unit(0.75, "cm")))
        
        plot = Heatmap(t(to_plot), 
        cluster_rows = F, 
        cluster_columns = F,
        show_row_names = F,
        show_column_names = T,
        row_title = "Cells",
        col = colorRamp2(seq(-5, 10, length.out=299), rev(colorRampPalette(brewer.pal(11, "RdBu"))(299))))
        pdf(paste0(plot_path, "NMF cluster ", cluster_name, ".pdf"), height = 12, width = 7.5)
        
        print(plot + ha)
        dev.off()
        
        return(to_plot)
}

#The function takes as argument a list of meta-signatures (a list of outputs from create.metasignature). It outputs a data frame of GO Terms that can help
#identify the cell types specific to each meta-signature.

metasignature.go.terms = function(metasignatures, results_path, num_terms = 20)
{ 
        dbs <- c("GO_Biological_Process_2015")
        go_terms = list()
        for (i in 1:length(metasignatures))
        { 
                 go_terms[[i]] = enrichr(rownames(metasignatures[[i]]),dbs)
        }      

        go_terms = data.frame(lapply(go_terms, function(x) x[[1]][1:num_terms,1]))
        colnames(go_terms) = paste0("MetaSignature", seq(1,length(metasignatures)))
        write.csv(go_terms, paste0(results_path, "Go Terms.csv"))
        return(go_terms)
}





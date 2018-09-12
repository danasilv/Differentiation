setwd("/Volumes/ahg_regevdata2/projects/Glioma_differentiation")

library(Seurat)
library(dplyr)
library(Matrix)
library(plotrix)
CITESEQ_COUNTS_THRESHOLD = 200

#Define sample to be analyzed, and the sub-folders for tpm, figures, and results.
sample = "MGG18"

tpm_path = paste0("TPM/", sample)
figures_path = paste0("figures/", sample)
results_path = paste0("results/", sample)

# Read in cell hashtags.
hashtag_files = list.files(tpm_path)
hashtag_files = hashtag_files[grep("cellnames", hashtag_files)]

hashtag = list()
for (i in 1:length(hashtag_files))
{hashtag[[i]] = read.table(paste0(tpm_path, "/", hashtag_files[i]), sep = ",", stringsAsFactors = F)}
hashtag = do.call(rbind, hashtag)
colnames(hashtag) = c("Barcode", "TimePoint")

#Read in tumour data.
tumor_10x = Read10X(paste0(tpm_path, "/GRCh38"))
tumor_data = CreateSeuratObject(raw.data = tumor_10x, min.cells = 3, min.genes = 200, project = sample)

#Generate QC plots.
mito_genes = grep(pattern = "^MT-", x = rownames(x = tumor_data@data), value = TRUE)
percent_mito = Matrix::colSums(tumor_data@raw.data[mito_genes,
                                                   ])/Matrix::colSums(tumor_data@raw.data)
tumor_data = AddMetaData(object = tumor_data, metadata = percent_mito, col.name
                         = "percent_mito")
pdf(file = paste0(figures_path,"/QC.pdf"), height = 6, width = 6)
VlnPlot(object = tumor_data, features.plot = c("nGene", "nUMI",
                                               "percent_mito"), nCol = 3, size.title.use = 14)
dev.off()

#Cell filtering.
tumor_data <- FilterCells(object = tumor_data, subset.names = c("nGene", "percent_mito"),
                          low.thresholds = c(2000, -Inf), high.thresholds = c(6000, 0.2))

#Add in hashtags as metadata.
hashtag_metadata = data.frame(Barcode = tumor_data@cell.names, stringsAsFactors = F)
hashtag_metadata = left_join(hashtag_metadata, hashtag, by = "Barcode")
hashtag_metadata$TimePoint[is.na(hashtag_metadata$TimePoint)] = "None"

timepoints = unique(hashtag_metadata$TimePoint)

for (i in 1:length(timepoints))
{
  temp = hashtag_metadata$TimePoint == timepoints[i]
  hashtag_metadata = cbind(hashtag_metadata, temp)
  colnames(hashtag_metadata)[ncol(hashtag_metadata)] = timepoints[i]
}
rownames(hashtag_metadata) = hashtag_metadata$Barcode
hashtag_metadata = dplyr::select(hashtag_metadata, -Barcode, -TimePoint)
hashtag_metadata[,1:ncol(hashtag_metadata)] = apply(hashtag_metadata[,1:ncol(hashtag_metadata)], 1:2, as.numeric)
tumor_data = AddMetaData(tumor_data, hashtag_metadata)

#Data normalization
tumor_data <- NormalizeData(object = tumor_data, normalization.method = "LogNormalize",
                        scale.factor = 10000)

##Find variable genes
pdf(file = paste0(figures_path,"/Variable Genes.pdf"), height = 6, width = 6)
tumor_data <- FindVariableGenes(object = tumor_data, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()

#Scale Data.
tumor_data <- ScaleData(tumor_data, do.center = TRUE, do.scale = FALSE, vars.to.regress = c("percent_mito"))

#Perform PCA.
tumor_data <- RunPCA(object = tumor_data, pc.genes = tumor_data@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 15)

#PCA plot.
pdf(paste0(figures_path, "/Principal Components.pdf"), width = 10, height = 12)
PCHeatmap(object = tumor_data, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()
 
#JackStraw
tumor_data <- JackStraw(object = tumor_data, num.replicate = 100, display.progress = FALSE)
pdf(paste0(figures_path, "/JackStraw plot.pdf"), width = 10, height = 12)
JackStrawPlot(object = tumor_data, PCs = 1:20)
dev.off()

#Find clusters.
tumor_data <- FindClusters(tumor_data, reduction.type = "pca", dims.use = 1:16,
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)

#TSNE.
tumor_data <- RunTSNE(object = tumor_data, dims.use = 1:16, do.fast = TRUE)

pdf(paste0(figures_path, "/TSNE Plot.pdf"), width = 6, height = 6)
TSNEPlot(object = tumor_data)
dev.off()

#look at time points
pdf(paste0(figures_path, "/Time Points Plot.pdf"), width = 10, height = 10)
FeaturePlot(object = tumor_data, features.plot = c("None", "CONTROL", "1H", "3H",
                                              "6H", "9H", "12H", "24H", "48H"),
            cols.use = c("#BEBEBE4D", "blue"),
            reduction.use = "tsne")
dev.off()

#look at expression of MHC genes.
pdf(paste0(figures_path, "/MHC1 genes Plot.pdf"), width = 8, height = 8)
FeaturePlot(object = tumor_data, features.plot = c("B2M", "HLA-A", "HLA-B", "HLA-C"),
            cols.use = c("#BEBEBE4D", "blue"),
            reduction.use = "tsne")
dev.off()

#Differential expression analysis.
markers <- FindAllMarkers(tumor_data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers, paste0(results_path, "/All Marker Genes.csv"))

#find GO terms (using EnrichR) of each cluster.
genes = list()
clusters = unique(markers$cluster)

for (i in 1:length(clusters))
{
     genes[[i]] = markers[ markers$cluster == clusters[i],]
}
go_terms = metasignature.go.terms(genes, num_terms = 20)
write.csv(go_terms, paste0(results_path, "/GO terms.csv"))

#GBM2.0 signatures.
gene_signatures = new.sigs (tumor_data)
tumor_data = AddMetaData(tumor_data, gene_signatures)

pdf(paste0(figures_path, "/Gene Signatures Plot.pdf"), width = 16, height = 12)
FeaturePlot(object = tumor_data, features.plot = colnames(gene_signatures),
            cols.use = c("#BEBEBE4D", "blue"),
            reduction.use = "tsne")
dev.off()

#Look at expression of some marker genes.
pdf(paste0(figures_path, "/Differentiation Genes Expression Plot.pdf"), width = 10, height = 8)
FeaturePlot(object = tumor_data, features.plot = c("SOX2", "NES", "HES1", "EGFR",
                                                   "GFAP", "RBFOX3", "TUBB3"),
            cols.use = c("#BEBEBE4D", "blue"),
            reduction.use = "tsne")
dev.off()

#Force timepoint into the "ident" section for differential gene expression analysis between timepoints.
hashtag_metadata = data.frame(Barcode = tumor_data@cell.names, stringsAsFactors = F)
hashtag_metadata = left_join(hashtag_metadata, hashtag, by = "Barcode")
hashtag_metadata$TimePoint[is.na(hashtag_metadata$TimePoint)] = "None"
rownames(hashtag_metadata) = hashtag_metadata$Barcode

tumor_data_timepoints = tumor_data
tumor_data_timepoints = AddMetaData(tumor_data_timepoints, hashtag_metadata)
tumor_data_timepoints = SetAllIdent(tumor_data_timepoints, "TimePoint")



#Compare all time points to control.
timepoints_to_compare = c("1H", "3H", "6H", "9H", "12H", "24H", "48H")
timepoint_comparisons = list()
for (i in 1:length(timepoints_to_compare))
{ 
        temp = FindMarkers(tumor_data_timepoints, timepoints_to_compare[i], "CONTROL", only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
        temp = data.frame(temp, Comparison = paste(timepoints_to_compare[i], "vs CONTROL"))
        timepoint_comparisons[[i]] = temp
}
        timepoint_comparisons_out = do.call(rbind, timepoint_comparisons)
        write.csv(timepoint_comparisons_out, paste0(results_path, "/Timepoint comparisons UP.csv"))
        
        
for (i in 1:length(timepoints_to_compare))
{ 
        temp = FindMarkers(tumor_data_timepoints, "CONTROL", timepoints_to_compare[i], only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
        temp = data.frame(temp, Comparison = paste("CONTROL vs", timepoints_to_compare[i]))
        timepoint_comparisons[[i]] = temp
}
        timepoint_comparisons_out = do.call(rbind, timepoint_comparisons)
        write.csv(timepoint_comparisons_out, paste0(results_path, "/Timepoint comparisons DOWN.csv"))


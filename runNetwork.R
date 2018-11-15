source("SingleCellClustering.R")
source("CNVDetection.R")
source("Network.R")


#Define file paths
tumor = "MGH143"
path.start <- "/Volumes/ahg_regevdata2/projects/Glioma_differentiation/" 
data.path = paste0(path.start,"alignreads/SS2/",tumor,"/RNAseq_pipeline/")
plot.path = paste0(path.start, "/figures/pathways/")
results.path = paste0(path.start,"/results/pathways/")
resources.path = "/Volumes/ahg_regevdata2/projects/Glioma_differentiation/resources/"

#data_files contains the filenames for single-cell data.
#smartseq_data stores all single-cell data.
#sample_ident stores the sample identity for all the cells contained in nucseq_data.

#Read in data, convert to dataframe.
data.files = paste0(data.path, list.files(data.path)[grepl("MGH",list.files(data.path))])
smartseq.data = lapply(data.files, function(x) read.table(x, sep = "\t", header = T, row.names = 1))
smartseq.data = do.call(cbind, smartseq.data)
smartseq.data[is.na(smartseq.data)] = 0
smartseq.data = smartseq.data[colnames(smartseq.data) != "X2transcript_id"]

#Read in tables of normal cells from Itay.
# macrophage = read.table(paste0(data_path, "/macrophage.csv"), sep = "\t", header = T)
# oligo = read.table(paste0(data_path, "/oligo.csv"), sep = "\t", header = T)
# macrophage = gsub("-", "_", macrophage$Macrophage)
# oligo = gsub("-", "_", oligo$Oligodendrocytes)

#Retrieve sample identity.
sample_ident = substr(colnames(smartseq.data), 1, 6)

#For finding normal cells.
#colnames(smartseq_data) = gsub("\\.", "_", colnames(smartseq_data))

#############A. Process data and cluster########################

#Process data.
smartseq.data = tpm.process(smartseq.data, sample_ident, plot.path, nGene_cutoff_low = 2500, 
                            nGene_cutoff_high = 10000, Ea_cutoff = 4)

sample.ident = smartseq.data$sample_ident
smartseq.data.TPM = smartseq.data$TPM

# Create Seurat
srtTumor <- CreateSeuratObject(raw.data = smartseq.data.TPM, min.cells = 3, min.genes = 400, project = tumor)

# Mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = rownames(x = srtTumor@data), value = TRUE, ignore.case = T)
percent.mito <- Matrix::colSums(srtTumor@raw.data[mito.genes, ]) / Matrix::colSums(srtTumor@raw.data)
srtTumor <- AddMetaData(object = srtTumor, metadata = percent.mito, col.name = "percent.mito")

# House keeping genes (list from Itay Tirosh)
hkgenes <- read.table("/Volumes/ahg_regevdata2/projects/Glioma_scGenetics/resources/tirosh_house_keeping.txt", skip = 2)
hkgenes <- as.vector(hkgenes$V1)
hkgenes.found <- which(toupper(rownames(srtTumor@data)) %in% hkgenes)  # remove hkgenes that were not found
n.expressed.hkgenes <- Matrix::colSums(srtTumor@data[hkgenes.found, ] > 0)
srtTumor <- AddMetaData(object = srtTumor, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

VlnPlot(object = srtTumor, features.plot = c("nGene"))
ggsave(paste0(plot.path, 'nGene.png'), width = 10, height = 3)
#VlnPlot(object = srtTumor, features.plot = c("percent.mito"))
#ggsave(paste0(plot.path, 'percentMt.png'), width = 10, height = 3)
VlnPlot(object = srtTumor, features.plot = c("n.exp.hkgenes"))
ggsave(paste0(plot.path, 'nexphkgenes.png'), width = 10, height = 3)

GenePlot(object = srtTumor, gene1 = "nGene", gene2 = "n.exp.hkgenes")
ggsave(paste0(plot.path, 'nGenevhkgenes.png'), width = 5, height = 4)

srtTumor <- NormalizeData(object = srtTumor, normalization.method = "LogNormalize", scale.factor = 100000)
srtTumor@data[is.na(srtTumor@data)] <- 0 # why do i need to do this? Shpouldn't normalize data take care of it?


# #Define colours for plotting
# colours = brewer.pal(length(unique(sample_ident)), "Spectral")
# names(colours) = unique(sample_ident)
# colours = list(Sample = colours)
# 
# #Keep normal cells only.
# malignant_data = smartseq_data[!(colnames(smartseq_data) %in% c(macrophage, olig
# 
# #start with one sample only.
# sample_ident = substr(colnames(malignant_data), 1, 6)
# test = malignant_data[sample_ident == "MGH143"]

#start with gene signatures from Neftel et. al.
gbm2.sig.mean <- gbm2.sigs(srtTumor)
srtTumor <- AddMetaData(srtTumor, gbm2.sig.mean ,colnames(gbm2.sig.mean))



# Score pathways from reactome
#pathway.scores <- calcPathwayScores(list(t(smartseq.data.TPM)), paste0(resources.path, "REACTOME/"))
#or
pathway.scores <- read.csv("/Volumes/ahg_regevdata2/projects/Glioma_differentiation/results/pathways/MGH143_reactome_enrichment.csv", sep = ",", row.names = 1)


# Calculate correlation between the signatures and the pathways
all.pathway.ordered.rownames <- order(rownames(as.data.frame(pathway.scores)))
all.pathway.enrichment.ordered <- as.data.frame(pathway.scores)[all.pathway.ordered.rownames,]

# Remove pathways with sd == 0
sd.pathways <- sapply(all.pathway.enrichment.ordered, sd, na.rm = TRUE) != 0
all.pathway.enrichment.ordered = all.pathway.enrichment.ordered[,sd.pathways]

for (sig.name in colnames(gbm2.sig.mean)){
 
  curr.sig.enrichment <- FetchData(srtTumor,sig.name)
  
  sig.ordered.rownames <- order(row.names(curr.sig.enrichment))
  curr.sig.enrichment.ordered <- curr.sig.enrichment[sig.ordered.rownames,]
  
  
  if (dim(curr.sig.enrichment)[1] != dim(all.pathway.enrichment.ordered)[1]) {
    print ("Error in the number of cells for the signature vs the number of cells from the pathways")
  } else {
  
    # Calc correlation between current signature and all pathways
    pathways.cor <- cor(curr.sig.enrichment.ordered, all.pathway.enrichment.ordered, method = "spearman")
    pathways.cor <- as.data.frame(t(pathways.cor))
    pathways.cor <- pathways.cor[order(-pathways.cor$V1), , drop = FALSE] # sort according to correlation
    
    # print
    write.csv(pathways.cor,paste0(results.path,sig.name,"_spearmancorr.csv"))
  }
}


# geneset = read.csv(paste0(resources.path, "GBM_signatures.csv"), stringsAsFactors = F)
# geneset = as.list(geneset)
# geneset = lapply(geneset, function(x) x[!is.na(x)])
# 
# 
# 
# #network propagation 
# network_hippie = load_hippie(paste0(resources.path, "hippie_current.txt"))
# network_hippie = network_hippie[network_hippie$score >= 0.6,1:2]
# network_biogrid = load_biogrid(paste0(resources_path, "BIOGRID-ALL-3.4.162.tab2.txt"))
# network_trrust = load_trrust(paste0(resources_path, "trrust_rawdata.human.tsv"))
# network_database = rbind(network_hippie, network_trrust)
# 
# #Build network
# #using adjacency matrix as edge weights.
# #TODO. Figure out the best way to construct the network
# 
# graph_out = build.network(geneset[[i]], malignant_data, network_database,  
#                           filter_steps = c(FALSE, 
#                                                                                             FALSE,  
#                                                                                             TRUE,  
#                                                                                             TRUE,  
#                                                                                             TRUE,  
#                                                                                             TRUE, 
#                                                                                             TRUE,  
#                                                                                             TRUE))
# W = create.w(graph_out)
# 
# 
# for (i in 1:length(geneset))
# { 
#         #Replace geneset in with the current geneset.
#         graph_out$P_0 = rep(0, length(graph_out$P_0))
#         graph_out$P_0 [graph_out$gene_names %in% geneset[[i]]] = 1
#         
#         # Normalize the P_0 to sum == 1
# 
#         #propagation.
#         P = propagate(graph_out, W, alpha = 0.6, convergence_cutoff = 1e-6)
#         network_significance = network.significance (P, graph_out, W, pval_cutoff = 0.10,  n_permutations = 100)
#         #TODO. the calculation of network_significance is very slow
#         
#         #write out the results
#         write.csv(network_significance, paste0(results_path, names(geneset)[i], ".csv"))
#         
#         # ToDo Dana: filter according to a list, such as TFs or enzymes
#         
#         #plotting.
#         #plot.graph(graph_out, P, file1 = paste0(plot_path, names(geneset)[i], " graph initial correlations.pdf"), 
#                   # file2 = paste0(plot_path, names(geneset)[i], " graph final correlations.pdf"))
# }
# #TODO. end a loop here.
# 
# #Projection of cells
# # Creating a graph where each node is a gene and each edge is from the PPI, returning a matrix of priors, where each prior is the transcriptomic state of the cell
# projection = project.cells(test, graph_out) 
# 
# #pick two cells.
# projection_1 = projection 
# projection_1$P_0 = projection_1$P_0[,10]
# 
# projection_2 = projection
# projection_2$P_0 = projection_2$P_0[,125]
# 
# 
# #Propagate a single cell.
# P_1 = propagate(projection_1, W)
# testa = P_1[[1]]
# testb = P_1[[21]]
# plot(testa, testb)
# P_1 = cbind(names(testa)[order(testa, decreasing=TRUE)],
# colnames(testb)[order(testb, decreasing = TRUE)])
# #Propagate another cell.
# P_2 = propagate(projection_2, W)
# testa = P_2[[1]]
# testb = P_2[[21]]
# plot(testa, testb)
# P_2 = cbind(names(testa)[order(testa, decreasing=TRUE)],
#             colnames(testb)[order(testb, decreasing = TRUE)])
# 
# trial = cbind(P_1, P_2)
# 
# 
# network_significance = network.significance(P, graph_out, W)
# write.csv(network_significance, "network significance adjacency matrix.csv")
# 
# #using gene-gene correlations as edge weights.
# graph_out = build.network(signatures$Mes, malignant_data, network_database, use_adjacency_matrix = F)
# W = create.w(graph_out)
# P = propagate(graph_out, W)
# plot.graph(graph_out, P, file1 = "graph initial correlations.pdf", file2 = "graph final correlations.pdf")
# network_significance = network.significance(P, graph_out, W)
# write.csv(network_significance, "network significance correlations.csv")
# 
# #Add in Tfs.
# network_database2 = trrust()
# network_database = list(network_database, network_database2)
# network_database = combine.network(network_database)
# graph_out = build.network(signatures$Mes, malignant_data, network_database, use_adjacency_matrix = F)
# W = create.w(graph_out)
# P = propagate(graph_out, W)
# plot.graph(graph_out, P, file1 = "graph initial correlations with TF.pdf", file2 = "graph final correlations with TF.pdf")
# network_significance = network.significance(P, graph_out, W)
# write.csv(network_significance, "network significance correlations with TF.csv")
# 
# network_significance_tf = network.significance.annotate.tf(network_significance, trrust(), filter = T)
# write.csv(network_significance_tf, "network significance correlations with TF TFs only.csv")
# 
# #REACTOME pathway propagation.
# reactome = reactome()
# trrust_nw = trrust()
# network = rbind(reactome, trrust_nw)
# #remove NA.
# network = network[!is.na(network[,1]) & !is.na(network[,2]),]
# 
# graph_out = build.network(geneset = 0, tpm = data, network_database = network,
#                               use_adjacency_matrix = TRUE,
#                               filter_steps = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE))
# W = create.w(graph_out)
# graph_out = project.cells(data, graph_out)
# 
# #
# P = propagate(graph_out, W)
# 
# P_plot = list()
# P_plot[[1]] = P[[1]][20,]
# P_plot[[2]] = P[[1]][100,]
# plot.graph(graph_out, P_plot)
# 
# plot.graph = function(graph_out, P, file1 = "graph initial.pdf", file2 = "graph final.pdf")
#         
# 

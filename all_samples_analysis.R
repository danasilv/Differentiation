setwd("~/Suva Lab/Dana Project 1/06-27-2018 All Samples Analysis")
source("~/Suva Lab/Dana Project 1/Functions/Functions.R")

#Read in data
data = read.csv("~/Suva Lab/Dana Project 1/TPM/plate.csv", row.names = 1)
exome_cnv_files = list.files("~/Suva Lab/Dana Project 1/EXOME")
exome_cnv_files = exome_cnv_files[grep("segtab", exome_cnv_files)]

#get sample identity.
sample_ident = unlist(lapply(strsplit(colnames(data), "_"), function(x) x[1]))
sample_ident [ grep("RMGH", sample_ident)] = "Reference"

#ID non-malignant cells.
data = tpm.process(data, hk, sample_ident)

#get sample identity, after filtering.
sample_ident = unlist(lapply(strsplit(colnames(data), "_"), function(x) x[1]))
sample_ident [ grep("RMGH", sample_ident)] = "Reference"

clusters = tpm.cluster(data, sample_ident)
cluster_genes = tpm.expressed.genes(data, clusters)

#Based on cluster_genes,
#Cluster 1 = malignant cells.
#Cluster 2 = immune cells.
#Cluster 3 = oligo cells.

hc_clusters = cutree(hc, k = 3)
malignant = hc_clusters == 1
immune = hc_clusters == 2
oligo = hc_clusters == 3

#Generate tsne plot
data_tsne = Rtsne(t(data), dims = 2, perplexity = 50, check_duplicates = F)

celltype = rep("Malignant", length(malignant))
celltype[immune] = "Immune"
celltype[oligo] = "Oligodendrocyte"
Sample = substr(colnames(data), 1,6)

toplot = data.frame(TSNE_1 = data_tsne$Y[,1], TSNE_2 = data_tsne$Y[,2], Cell_type = celltype, Sample = Sample)

plot = ggplot(toplot, aes(x = TSNE_1, y = TSNE_2, colour = Cell_type)) + geom_point()
plot = plot + theme_classic()
plot = plot + labs (color = "Cell type")
plot = plot + scale_color_manual(values=c("#E69F00", "#999999",  "#56B4E9"))
pdf("TSNE by cell type.pdf", height = 6, width = 6)
print(plot)
dev.off()

plot = ggplot(toplot, aes(x = TSNE_1, y = TSNE_2, colour = Sample)) + geom_point()
plot = plot + theme_classic()
plot = plot + labs (color = "Sample")
plot = plot + scale_color_manual(values=c("#E69F00", "#FF00FF",  "#56B4E9", "#0000FF", "#00FFFF", "#FF0000"))
pdf("TSNE by Sample.pdf", height = 6, width = 6)
print(plot)
dev.off()



#Run NMF.
malignant_data = data[,malignant]
#Only 1 "Reference MGH105" cell got placed into the malignant cell cluster.
sample = substr(colnames(malignant_data),1,6)
unique_samples = unique(sample)
malignant_data = malignant_data[,sample != "RMGH10"]
nmf_signatures = tpm.to.nmf(malignant_data)
signature_scores = calc.signature.scores(malignant_data, nmf_signatures)
signature_cluster = cluster.signature.scores(signature_scores)

#Create meta-signature indices. 
k = 6
signature_indices = cutree(signature_cluster, k = k)
signature_indices_listed = list()

for (i in 1:k) signature_indices_listed[[i]] = signature_indices == i

#Create meta-signatures.
meta_signatures = list()
for (i in 1:k) meta_signatures[[i]] = create.metasignature(malignant_data, nmf_signatures, signature_scores, 
                                                           signature_indices_listed[[i]], paste0 ("MetaSignature ",i))
go_terms = metasignature.go.terms(meta_signatures, num_terms = 20)
 
#InferCNV.
Ecnv_smoothed = infer.cnv(data, gencode)
genes_cn = match.exome.cnv(Ecnv_smoothed, exome_cnv_files)
malignant_cn = id.malignant(Ecnv_smoothed, genes_cn)

#compare malignant cells defined by CNV, and malignant cells defined by clustering.
malignant_exprs = colnames(data)[malignant]
length(intersect(malignant_cn, malignant_exprs))
#For now, just use malignant cells as defined by gene expression.
oligo = colnames(data)[oligo]

Ecnv_smoothed_ref = baseline.correction(Ecnv_smoothed, malignant_exprs, oligo, method = "oligo", output = "malignant")
cnv_plot = plot.cnv(Ecnv_smoothed_ref, genes_cn)
cnv_plot   
        
#network propagation 
network_database = string()

#using adjacency matrix as edge weights.
graph_out = build.network(signatures$Mes, malignant_data, network_database)
W = create.w(graph_out)
P = propagate(graph_out, W)
plot.graph(graph_out, P, file1 = "graph initial adjacency matrix.pdf", file2 = "graph final adjacency matrix.pdf")
network_significance = network.significance(P, graph_out, W)
write.csv(network_significance, "network significance adjacency matrix.csv")

#using gene-gene correlations as edge weights.
graph_out = build.network(signatures$Mes, malignant_data, network_database, use_adjacency_matrix = F)
W = create.w(graph_out)
P = propagate(graph_out, W)
plot.graph(graph_out, P, file1 = "graph initial correlations.pdf", file2 = "graph final correlations.pdf")
network_significance = network.significance(P, graph_out, W)
write.csv(network_significance, "network significance correlations.csv")

#Add in Tfs.
network_database2 = trrust()
network_database = list(network_database, network_database2)
network_database = combine.network(network_database)
graph_out = build.network(signatures$Mes, malignant_data, network_database, use_adjacency_matrix = F)
W = create.w(graph_out)
P = propagate(graph_out, W)
plot.graph(graph_out, P, file1 = "graph initial correlations with TF.pdf", file2 = "graph final correlations with TF.pdf")
network_significance = network.significance(P, graph_out, W)
write.csv(network_significance, "network significance correlations with TF.csv")

network_significance_tf = network.significance.annotate.tf(network_significance, trrust(), filter = T)
write.csv(network_significance_tf, "network significance correlations with TF TFs only.csv")

#REACTOME pathway propagation.
reactome = reactome()
trrust_nw = trrust()
network = rbind(reactome, trrust_nw)
#remove NA.
network = network[!is.na(network[,1]) & !is.na(network[,2]),]

graph_out = build.network(geneset = 0, tpm = data, network_database = network,
                              use_adjacency_matrix = TRUE,
                              filter_steps = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE))
W = create.w(graph_out)
graph_out = project.cells(data, graph_out)

#
P = propagate(graph_out, W)

P_plot = list()
P_plot[[1]] = P[[1]][20,]
P_plot[[2]] = P[[1]][100,]
plot.graph(graph_out, P_plot)

plot.graph = function(graph_out, P, file1 = "graph initial.pdf", file2 = "graph final.pdf")
        


#####################Load in R libraries#########################
library(Seurat)
library(Matrix)
library(stringr)
library(DiagrammeR)
library(fifer)
library(NMF)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamps)
library(circlize)
library(ggplot2)
library(enrichR)
library(cluster)
library(tidyr)
library(mygene)
library(dplyr)
library(igraph)
library(Rtsne)
library(scales)
library(biomaRt)
library(diptest)

#An example of how to use these functions:
#network_database2 = trrust()
#network_database = list(network_database, network_database2)
#network_database = combine.network(network_database)
#graph_out = build.network(signatures$Mes, malignant_data, network_database, use_adjacency_matrix = F)
#W = create.w(graph_out)
#P = propagate(graph_out, W)
#plot.graph(graph_out, P, file1 = "graph initial correlations with TF.pdf", file2 = "graph final correlations with TF.pdf")
#network_significance = network.significance(P, graph_out, W)
#write.csv(network_significance, "network significance correlations with TF.csv")

#network_significance_tf = network.significance.annotate.tf(network_significance, trrust(), filter = T)
#write.csv(network_significance_tf, "network significance correlations with TF TFs only.csv")


#####################Load in dependencies#########################

hk = read.table("Resources/tirosh_house_keeping.txt", skip = 2)
signatures = read.table("Resources/GBM_signatures.csv", header = TRUE, sep = ",", stringsAsFactors = F)
signatures = as.list(signatures)
signatures = lapply(signatures, function(x) x[!is.na(x)])
gencode = read.table("Resources/gencode_v19_gene_pos.txt")

###################Network propagation############################
#TODO: for all of these, add a column which specifies whether the edge should be one direction or both directions.

#This function loads in the biogrid database for network propagation. To use, network_database = load_biogrid(path to database file).
load_biogrid = function(biogrid_path = "~/Suva Lab/Dana project 2/RESOURCES/BIOGRID-All-3.4.162.tab2.txt")
{
        #read in database.
        network_database = read.delim(biogrid_path,
                                      header = T, row.names = 1, stringsAsFactors = F)
        #keep interactions for human.
        network_database = network_database[network_database$Organism.Interactor.A == 9606 & network_database$Organism.Interactor.A == 9606,]
        network_database = dplyr::select(network_database, Official.Symbol.Interactor.A, Official.Symbol.Interactor.B)
        colnames(network_database) = c("from", "to")
        #remove cases where a gene "interacts with itself".
        network_database = network_database[network_database$from != network_database$to,]
        return (network_database)
}

#This function loads in the trrust database of TF-target gene interactions for network propagation. To use, network_database =
#load_trrust(path to database file)
load_trrust = function(trrust_path = "~/Suva Lab/Dana Project 2/RESOURCES/trrust_rawdata.human.tsv" )
{
        network_database = read.table(trrust_path, sep = "\t")
        colnames(network_database) = c("from", "to", "type", "ID")
        network_database = dplyr::select(network_database, from, to)
}

#This function loads in the hippie database of protein-protein interactions for network propagation. To use, network_database = 
#load_hippie(path to database file)
load_hippie = function(hippie_path = "~/Suva Lab/Dana Project 2/RESOURCES/hippie_current.txt")
{
        network_database = read.delim(hippie_path, header = F, stringsAsFactors = F)
        network_database = select(network_database, V1, V3, V5)
        colnames(network_database) = c("from", "to", "score")
        #remove self-interactions.
        network_database = network_database[network_database$from != network_database$to,]
        network_database[,1:2] = apply(network_database[,1:2], 1:2, function(x) gsub("_HUMAN", "", x))
        return(network_database)
}

#This function loads in the STRING database of protein-protein interactions for network propagation. To use, network_database = 
#load_string(path to database file)
load_string = function(string_path = "~/Suva Lab/Dana Project 2/RESOURCES/9606.protein.links.full.v10.5.txt", score_cutoff = 400)
{
        network_database = read.delim(string_path, header = T, stringsAsFactors = F, sep = " ")
        
        #Filtering and re-formating. 
        network_database = network_database[network_database$combined_score >= score_cutoff,]
        network_database = dplyr::select(network_database, protein1, protein2, combined_score)
        colnames(network_database) = c("from", "to", "score")
        
        #convert the STRINGdb IDs to gene IDs.
        mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
        nodes = gsub("9606.", "", unique(c(network_database$from, network_database$to)))
        query = getBM(mart = mart,attributes=c("ensembl_peptide_id","hgnc_symbol","entrezgene"),filters = "ensembl_peptide_id",values = nodes)
        
        query$ensembl_peptide_id = paste0("9606.", query$ensembl_peptide_id)
        
        #Is there a faster way to find and replace? the method below is faster than iterating over "query" and replacing from "network_database".
        #from.
        test = split.data.frame(network_database, network_database$from)
        from_replace_ID = sapply(test, function(x) unique(x$from))
        from_replace_gene = as.character(sapply(from_replace_ID, function(x) query$hgnc_symbol[query$ensembl_peptide_id == x]))
        
        n = sapply(test, function(x) nrow(x))
        from_replace = vector()
        for (i in 1:length(n)) from_replace = c(from_replace, rep(from_replace_gene[i], n[i]))
        
        test = do.call('rbind', test) 
        test$from = from_replace
        network_database = test
        
        #to.
        test = split.data.frame(network_database, network_database$to)
        to_replace_ID = sapply(test, function(x) unique(x$to))
        to_replace_gene = as.character(sapply(to_replace_ID, function(x) query$hgnc_symbol[query$ensembl_peptide_id == x]))
        
        n = sapply(test, function(x) nrow(x))
        to_replace = vector()
        for (i in 1:length(n)) to_replace = c(to_replace, rep(to_replace_gene[i], n[i]))
        
        test = do.call('rbind', test) 
        test$to = to_replace
        network_database = test
        
        rownames(network_database) = 1:nrow(network_database)
        
        #remove unmatched.
        network_database = network_database[network_database$from %in% query$hgnc_symbol & network_database$to %in% query$hgnc_symbol,]
        
        #the score is not really a measure of interaction strength, rather a measure of confidence in the interaction.
        network_database = dplyr::select(network_database, -score)
        return(network_database)
}

#This function loads in the REACTOME database for network propagation. To use, network_database = load_reactome()
load_reactome = function(reactome_path = "~/Suva Lab/Dana Project 2/RESOURCES/REACTOME")
{
        #in mappings, second column is entrez gene ID, third column is gene symbol.
        mappings = read.delim(paste0(reactome_path, "/Homo_sapiens.gene_info"), header = F,
                              stringsAsFactors = F, sep = "\t", skip = 1)[,1:3]
        colnames(mappings) = c("Species", "Entrez", "Symbol")
        mappings = dplyr::filter(mappings, Species == 9606)
        
        
        interactions = read.delim(paste0(reactome_path, "/reactome.homo_sapiens.interactions.tab-delimited.txt"), header = T,
                                  stringsAsFactors = F)
        
        
        interactions = interactions[,c(3,6)]
        colnames(interactions) = c("from", "to")
        interactions = apply(interactions, 1:2, function(x) gsub("entrezgene/locuslink:", "", x))
        #Anything that isn't mapped to a gene symbol should be removed.
        interactions = apply(interactions, 1:2, function(x) as.numeric(x))
        interactions = interactions[!is.na(interactions[,1]) & !is.na(interactions[,2]), ]
        
        #convert Entrez to GENE SYMBOL.
        interactions = apply(interactions, 1:2, function(x) mappings$Symbol[mappings$Entrez == x][1])
        
        #remove self-interactions.
        interactions = interactions[interactions[,1] != interactions[,2],]
        interactions = data.frame(interactions, stringsAsFactors = F)
        
        return(interactions)
}

#This function takes as input a list of networks (from e.g. biogrid(), trrust()...), and combines them into one network.
#Note that duplicate removal will occur in the build.network() function.
combine.network = function(networks)
{
        networks = lapply(networks, function(x) x[,1:2])
        networks = do.call("rbind", networks)
        return(networks)
}

#This function is called within the build.network function to check connectivity of the network and remove islands.
#TODO: for fully connected networks, have to be able to go from *any* gene to all other genes.

fix.connectivity = function(network_keep)
{ 
        # network_keep[,1] = as.character(network_keep[,1])
        # network_keep[,2] = as.character(network_keep[,2])
        unique_genes = unique(c(network_keep[,1], network_keep[,2]))
        continue = TRUE
        
        while(continue)
        { 
                #Sample a gene.
                temp_gene = sample(unique_genes, 1)
                temp_network = temp_gene 
                length = vector()
                length[1] = 1
                
                #Get all connections to the sampled gene.
                for (i in 2:length(unique_genes))
                {         
                        temp_interactions = network_keep[network_keep$from %in% temp_network | network_keep$to %in% temp_network,]
                        temp_network = c(temp_gene,temp_interactions$from, temp_interactions$to)
                        temp_network = unique(temp_network)
                        length[i] = length(temp_network)
                        if (length[i] == length[i-1]) break
                }
                
                #Remove islands, if there are any.
                if (length(temp_network) == length(unique_genes))
                {
                        print("Network graph is already connected.")
                        continue = FALSE
                }else{
                        print("Network graph is not connected. Removing islands.")
                        #If the length of the current network is less than half, it is in the minority 
                        #and should be removed. If the length of the current network is more than half,
                        #it is in the majority and should be kept. 
                        if (length(temp_network) < length(unique_genes)/2)
                        {
                                genes_keep = unique_genes[!unique_genes %in% temp_network]   
                        }else{
                                genes_keep = unique_genes[unique_genes %in% temp_network]        
                        }
                        #remove unconnected genes.
                        network_keep = network_keep[network_keep$from %in% genes_keep | network_keep$to %in% genes_keep, ]  
                        #reset the unique genes list.
                        unique_genes = unique(c(network_keep[,1], network_keep[,2]))
                        
                        #If the network now contains all genes, can stop.
                        if (length(temp_network) == length(unique_genes))
                        {
                                print("Network graph is now connected.")
                                continue = FALSE
                        }
                }
        }
        return(network_keep)
}

#TODO: self loops?
#TODO: for one-directional networks, how to combine with bi-directional networks. especially for step 5, it should only be
#Applied to the bidirectional part.

#This function takes as arguments a geneset (list of genes), filtered tpm, and a network database (e.g. from biogrid()). 
#It will return a list object which includes all the genes within the supplied geneset and their genetic interactions from biogrid.
#The network will only contain genes that are actually detected within the single-cell rna-seq data, and interactions of gene pairs
#only if the two genes are positively correlated within the tpm dataset. The nodes will be annotated with "start", which = 1 for
#genes in the geneset. The node edges will be annotated with "weight", which is the either the gene pair correlation values normalized to 1
#(if use_adjacency_matrix = FALSE), or just 1/degree of the node (if use_adjacency_matrix = TRUE).
build.network = function(geneset, tpm, network_database, correlation_threshold = 0.1,
                         use_adjacency_matrix = TRUE,
                         filter_steps = c(TRUE, #Step 1. keep only interactions with genes in geneset.
                                          TRUE, #Step 2. keep only genes detected in tpm.
                                          TRUE, #Step 3. include all interactions within network.
                                          TRUE, #Step 4. remove duplicates
                                          TRUE, #Step 5. add back B -> A for A -> B. (do not do this step for directed networks)
                                          TRUE, #Step 6. remove duplicates again.
                                          TRUE, #Step 7. filter by correlation values
                                          TRUE))#Step 8. check network connectivity. (do not do this step for trrust only)
{
        network_keep = network_database
        #Step 1. search for genetic interactions which include the genes within the geneset.
        if (filter_steps[1] == TRUE)
        { 
                network_keep = network_database[network_database$from %in% geneset | network_database$to %in% geneset,]
        }
        #Step 2. keep only genes which are detected within the tpm.
        #alternative is to not do this filtering step, but use edge weights from another source.
        if (filter_steps[2] == TRUE)
        { 
                tpm_genes = rownames(tpm)
                network_keep = network_keep[network_keep$from %in% tpm_genes & network_keep$to %in% tpm_genes, ]
                rm(tpm_genes)
        }
        
        #Step 3. add in all other interactions between genes that are now included in the set.
        if (filter_steps[3] == TRUE)
        { 
                unique_genes = unique(c(network_keep[,1], network_keep[,2]))
                network_keep = rbind(network_keep, 
                                     network_database[network_database$from %in% unique_genes & network_database$to %in% unique_genes,])
                rm(unique_genes)
        }
        
        #Step 4. remove duplicates, and cases where both A -> B and B -> A are present.
        #Because edge weights in two directions is not the same after normalization, every B -> A is added back again.
        if (filter_steps[4] == TRUE)
        { 
                check_dups = vector()
                for (i in 1:nrow(network_keep))
                {
                        #the sorting here converts B -> A to A -> B for removal of duplicates.
                        temp = sort(as.character(unlist(network_keep[i,1:2])))
                        temp = paste(temp[1],temp[2], sep =" ")
                        check_dups[i] = temp
                }
                network_keep = network_keep[!duplicated(check_dups),]
                rm(check_dups)
        }
        
        #Step 5.  add back the B -> A cases for each of the A -> B. This is because the weights in the two directions
        #will be different, after normalization.
        if (filter_steps[5] == TRUE)
        { 
                temp = data.frame(from = network_keep$to, to = network_keep$from)
                network_keep = rbind(network_keep, temp)
                rm(temp)
        }
        
        #TODO: why does this happen??
        
        #Step 6. remove duplicates again, because the previous step would have introduced new duplicates.
        if (filter_steps[6] == TRUE)
        { 
                check_dups = vector()
                for (i in 1:nrow(network_keep))
                {
                        temp = as.character(network_keep[i,1:2])
                        temp = paste(temp[1],temp[2], sep ="")
                        check_dups[i] = temp
                }
                network_keep = network_keep[!duplicated(check_dups),]
                rm(check_dups)
        }
        
        #Determiniation of edge weights. The first case here is if we just want to use an adjacency matrix as edge weights (0,1), normalized by the degree.
        if (use_adjacency_matrix == TRUE)
        {
                print("Using adjacency matrix as edge weights, normalized by node degree.")
                network_keep = data.frame(network_keep, correlation = 1)
                #If scores are not provided. Generate gene correlations from tpm, use correlation values in edge weights.        
        }else if (ncol(network_keep) == 2){ 
                print("Interaction scores are not provided, building correlation matrix.")
                
                cor_matrix = cor(t(tpm), method = "pearson")
                network_keep = data.frame(network_keep, correlation = 0)
                
                for (i in 1:nrow(network_keep))
                {
                        network_keep$correlation[i] = cor_matrix[rownames(cor_matrix) == network_keep$from[i], 
                                                                 colnames(cor_matrix) == network_keep$to[i]]
                }
        }else{
                #TODO: needs to be fixed because adding a third column for edge direction.
                print("Using pre-supplied interaction scores as edge weights.")
                colnames(network_keep)[3] = "correlation"        
        }
        
        #Step 7. Keep only correlations above a specific threshold.
        if (filter_steps[7] == TRUE)
        {network_keep = network_keep[network_keep$correlation > correlation_threshold,]}
        
        #Step 8. Check connectivity of graph. If not connected, remove islands (under construction).
        if (filter_steps[8] == TRUE)
        { 
                network_keep = fix.connectivity (network_keep)
        }
        
        #Conversion of gene expression correlations to edge weights (decide on how to create weights)
        #All edge weights for a given node are divided by the degree of the node.
        #NOTE: this is now done within the "propagate" function.
        #network_keep = split(network_keep, network_keep$to)
        #network_keep = lapply(network_keep, function(x) data.frame(x, weight = x$correlation/nrow(x)))
        #network_keep = do.call(rbind, network_keep)
        #rownames(network_keep) = 1:nrow(network_keep)
        network_keep = group_by(network_keep, from) %>%
                mutate(weight = correlation/sum(correlation)) %>% ungroup()
        network_keep = data.frame(network_keep)
        
        #generate nodes
        genes = unique(unlist(network_keep[,1:2]))
        genes = data.frame(gene = genes, start = genes%in% geneset)
        
        #remove NAs.
        #returning the graph as a list of data frames.
        graph_out = list()
        #first object is the gene names.
        graph_out [[1]] = genes$gene[!is.na(genes$gene)]
        #second object is the propagation start sites.
        graph_out [[2]] = genes$start
        #third object is the edge weights.
        graph_out [[3]] = network_keep[!is.na(network_keep$from) & !is.na(network_keep$to),]
        names(graph_out) = c("gene_names", "P_0", "edge_weights")
        
        return(graph_out)
}

#This function converts the edge weights into a matrix format for network propagation.
create.w = function(graph_out)
{
        edge_weights = graph_out$edge_weights
        gene_names = graph_out$gene_names
        edge_weights[,1:2] = apply(edge_weights[,1:2], 1:2, function(x) which(gene_names == x))
        
        #convert the edge_weights into matrix form.
        W = data.matrix(matrix(0, nrow = length(gene_names), ncol = length(gene_names)))
        rownames(W) = colnames(W) = gene_names
        
        #Adjacency matrix.
        for (i in 1:nrow(edge_weights))
        {
                print(i)
                W[edge_weights$from[i], edge_weights$to[i]] = edge_weights$weight[i]
        }
        
        return (W)      
}

#Given a filtered tpm and the output of build.network, this function will output a matrix object containing the P_0 for 
#every cell. In the output, the rows are cells, and the columns are genes.
project.cells = function(tpm, graph_out)
{
        #TPM must be transformed such that it only contains positive values,
        #normalize across cell (columns)
        tpm_norm = apply(tpm, 2, function(x) (x-min(x))/(max(x)-min(x)))
        
        
        #Retrieve gene expression per cell.
        P_0_list = list()
        for (i in 1:ncol(tpm_norm))
        {
                print(i)
                temp_tpm = data.frame(Gene = rownames(tpm_norm), tpm_norm[,i])
                temp_P_0 = data.frame(Gene = graph_out$gene_names, Signal = rep(0, length(graph_out$gene_names)))
                for (j in 1:nrow(temp_P_0)) 
                {
                        if (as.character(temp_P_0[j,1]) %in% temp_tpm$Gene)
                        { 
                                temp_P_0[j,2] = temp_tpm[,2][temp_tpm$Gene == as.character(temp_P_0[j,1])]
                        }
                }
                temp_P_0 = as.vector(temp_P_0[,2])
                P_0_list[[i]] = temp_P_0
        }
        
        
        #Convert to matrix format, for matrix multiplication.
        #Rows are cells, columns are signals.
        P_0 = do.call(rbind, P_0_list)
        rownames(P_0) = colnames(tpm)
        colnames(P_0) = graph_out$gene_names
        
        #Normalize so that sum across each cell is 1.
        P_0 = apply(P_0, 1, function(x) x/sum(x))
        
        graph_out$P_0 = P_0
        return(graph_out)
}

#This function takes as argument a list generated by build.network(), and the edge weight matrix generated by create.w().
#It will return a list object that contains the values after each iteration of network propagation, with P[[1]] being
#the starting signal and P[[n]] being the result at convergence. It's expected that the P_0 which is supplied is already
#normalized to 1.

propagate = function(graph_out, W, n_iter = 1000, convergence_cutoff = 1e-12, alpha = 0.3, verbose = TRUE)
        #propagate with alpha = 0, 1 for p-val calculation.
        #TODO: Add a check for normalization.
        #TODO: normally would take 100s of iterations - is the network too connected?
        #TODO: how many genes? how many edges?
{ 
        #Inform status.
        if (verbose) print(paste0("Running network propagation for ", n_iter, 
                                  " iterations, or until convergence at a cutoff of ", convergence_cutoff,
                                  "."))
        #P_0 will be a logical vector for the geneset propagation. P_0 will be a numerical matrix
        #for the cell propagation.
        P = list()
        P_0 = graph_out$P_0
        if (is.logical(P_0)) {P[[1]] = as.numeric(P_0)
        vector = TRUE} else{
                P[[1]] = P_0
                vector = FALSE
        } 
        
        #Propagation. 
        #TODO. For cell propagation, each cell needs to be done separately. Otherwise some cells
        #will converge well before others do!
        for (i in 1:n_iter)
        {
                print(i)
                P[[i+1]] = alpha*P[[i]] + (1 - alpha) * (P[[i]] %*% W)
                #TODO: fix this - not needed.
                # if (vector) P[[i]] = as.vector(P[[i+1]])
                
                #TODO: normally would not use length.
                if (sum((P[[i+1]] - P[[i]])^2)/length(P[[i+1]]) <= convergence_cutoff) break
        }
        
        #TODO: results need to be normalized after propagation, divided by sum per cell? or is this necessary?
        if (verbose) print(paste0("Completed ", length(P) - 1, " iterations."))
        return (P)
}

#This function determines the significance of each gene identified by network propagation using a permutation test.
#It will return a data frame containing gene names, final values after propagation, fold difference between observed
#and permutation, and p-values for genes that have a p value less than the specified p value cutoff.

network.significance = function(P, graph_out, W, pval_cutoff = 0.05,  n_permutations = 1000)
{
        gene_names = graph_out$gene_names
        n_genes = sum(graph_out$P_0)
        
        #for storage of distribution
        observed = unlist(tail(P, 1))
        distributions = list()
        
        for (i in 1:n_permutations)
        {
                #select random genes.
                temp_genes = sample (gene_names, n_genes)
                temp_P_0 = gene_names %in% temp_genes
                temp_graph_out = graph_out
                temp_graph_out$P_0 = temp_P_0 
                
                #run propagation on these random genes.
                temp_P = tail(propagate(temp_graph_out, W, verbose = FALSE), 1)
                distributions[[i]] =  unlist(temp_P) 
        }
        
        distributions = as.data.frame(distributions)
        
        rownames(distributions) = gene_names
        colnames(distributions)[1:n_permutations] = paste("Permutation", seq(1, n_permutations))
        
        p = rowSums(apply(distributions, 2, function(x) x >= observed))/n_permutations
        foldchange = observed/rowMeans(distributions)
        
        network_significance = data.frame(observed = observed, foldchange = foldchange, p = p)
        rownames(network_significance) = gene_names
        network_significance = mutate(network_significance, gene = rownames(network_significance)) %>%
                arrange(desc(observed)) %>% 
                filter(p < pval_cutoff) %>%
                arrange(desc(foldchange))
        
        return(network_significance)
}  

#This function adds another column to the output of network.significance, indicating whether the gene is a regulator
#from the trrust database. "Filter" indicates whether only the regulator genes should be kept in the returned data frame.
network.significance.annotate.tf = function(network_significance, trrust, filter = TRUE)
{
        regulators = unique(trrust$from)
        network_significance = mutate(network_significance, regulator = gene %in% regulators)
        if (filter)
        {
                network_significance = filter(network_significance, regulator == TRUE)
                network_significance = select(network_significance, -regulator)
                return(network_significance)
        }else {return (network_significance)}
}

#This function makes a plot of a network, using the output of build.network (graph_out) and propagate (P).
plot.graph = function(graph_out, P, file1 = "graph initial.pdf", file2 = "graph final.pdf")
        
{
        graph = graph_from_data_frame(graph_out$edge_weights, directed = FALSE,vertices = as.character(graph_out$gene_names))
        layout = layout_as_tree(graph)
        
        colours = colorRampPalette(c('white', "darkgoldenrod4"))
        
        #Initial.
        cols = colours(100)[cut(P[[1]], breaks = 100)]
        pdf(file1, height = 8, width = 8)
        plot(graph, vertex.size = 4, vertex.color = adjustcolor(cols, alpha = 0.85),
             layout = layout, edge.width = graph_out$edge_weights$weight, vertex.label = NA)
        dev.off()
        
        #Final.
        cols = colours(100)[cut(unlist(tail(P,1)), breaks = 100)]
        pdf(file2, height = 8, width = 8)
        plot(graph, vertex.size = 4, vertex.color = adjustcolor(cols, alpha = 0.85),
             layout = layout, edge.width = graph_out$edge_weights$weight, vertex.label = NA)
        dev.off()
}

#This function takes the output of propagate, and calculates pathway scores using reactome, before
#and after network propagation.
#TODO: Clean it up.
pathway.scores = function(P, reactome_path = "~/Suva Lab/Dana Project 2/RESOURCES/REACTOME")
{
        pathways = read.delim(paste0(reactome_path, "/ReactomePathways.txt"), header = F,
                              stringsAsFactors = F)
        colnames(pathways) = c("Pathway", "Description", "Species")
        pathways = dplyr::filter(pathways, Species == "Homo sapiens")
        
        pathway_memberships = read.delim(paste0(reactome_path, "/Ensembl2Reactome_PE_All_Levels.txt"), header = F,
                                         stringsAsFactors = F)
        
        pathway_memberships = dplyr::filter(pathway_memberships, V8 == "Homo sapiens")
        pathway_memberships = pathway_memberships[,3:4]
        colnames(pathway_memberships) = c("Gene", "Pathway")
        
        pathway_memberships = dplyr::mutate(pathway_memberships, Symbol = "ABCDE")
        
        Genes = colnames(P[[1]])
        for (i in 1:length(Genes))
        {
                print(i)
                pathway_memberships$Symbol[grepl(Genes[i], pathway_memberships$Gene)] = Genes[i]
        }
        
        
        pathway_memberships = pathway_memberships[pathway_memberships$Symbol != "ABCDE",]
        
        #generate pathway scores for before and after propagation.
        before = P[[1]]
        after = P[[length(P)]]
        pathways = unique(pathway_memberships$Pathway)
        
        after_out = before_out = matrix(nrow = nrow (P[[1]]), ncol = length(pathways))
        rownames(after_out) = rownames(before_out) = rownames(P[[1]])
        colnames(after_out) = colnames(before_out) = pathways
        
        for (i in 1:ncol(after_out))
        {
                print(i)
                temp = pathway_memberships[pathway_memberships$Pathway == pathways[i],]
                temp = data.frame(Gene = unique(temp$Symbol))
                
                temp_P = t(P[[length(P)]])
                temp_P = data.frame(temp_P, Gene = rownames(temp_P))
                
                temp = inner_join(temp, temp_P, by = "Gene")
                temp = colMeans(temp[,2:ncol(temp)])
                
                after_out[,i] = temp
        }
        
        variance = apply(before_out, 2, var)
        test = before_out[,variance >= 0.05]
        
        for (i in 1:ncol(before_out))
        {
                print(i)
                temp = pathway_memberships[pathway_memberships$Pathway == pathways[i],]
                temp = data.frame(Gene = unique(temp$Symbol))
                
                temp_P = t(P[[1]])
                temp_P = data.frame(temp_P, Gene = rownames(temp_P))
                
                temp = inner_join(temp, temp_P, by = "Gene")
                temp = colMeans(temp[,2:ncol(temp)])
                
                before_out[,i] = temp
        }
        
        write.csv(before_out, "before_out.csv")
        write.csv(after_out, "after_out.csv")
}


####Functions for the analysis of single-cell RNA-Seq data 2.#####
#UPDATED: August 22, 2018. KL.

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


#####################Load in dependencies#########################

####################Phylogenetics############################


#This function is used in the parse.dot function below.
numextract = function(string){ 
        str_extract(string, "\\-*\\d+\\.*\\d*")
} 

#This function parses the .gv file format, and returns a dataframe containing all the paths and the last node for each path.
parse.dot = function (gvfile){
        
        #In the .gv file, toss the first line, and the last 3 lines. 
        gvfile = gvfile[2:(length(gvfile)-3)]
        gvfile = rev(gvfile)
        
        n_lines = length(gvfile)
        #Parse Edges.
        edges = gvfile[grep("->", gvfile)]
        edges = gsub("[^0-9 ]","",edges)
        temp_a = as.numeric(unlist(lapply(strsplit(edges, " "), function(x) x[1])))
        temp_b = as.numeric(unlist(lapply(strsplit(edges, " "), function(x) x[3])))
        edges = data.frame(From = temp_a, To = temp_b)
        
        #Parse labels.
        labels = gvfile[grep("label", gvfile)]
        parsed_labels = vector(mode = "character", length = length(labels))
        numbers = as.numeric(numextract(labels))
        genes = unlist(lapply(strsplit(labels, "label=\""), function(x) x[2]))
        genes = gsub("[^A-Za-z0-9 ]","",genes)
        genes = gsub(" style  filled","",genes)
        
        for (i in 1:length(numbers))
        {parsed_labels[numbers[i]+1] = genes[i]}
        
        #Replace the numbers in the edges data frame with labels.        
        edges = apply(edges, 1:2, function(x) parsed_labels[x+1])
        
        #Create paths.
        paths = character()
        
        for (i in 1:nrow(edges))
        {
                #for each entry in the paths collection, check the terminal node against the "From" column of the current edge.
                last_node = unlist(lapply(strsplit(paths, ":"), function(x) x[length(x)]))
                check = grep(edges[i,1], last_node)
                
                #If the "From" of the edge is not currently present, create a new entry in the paths list.
                if (length(check) == 0)
                {
                        paths = c(paths, paste0(edges[i,1],":",edges[i,2]))
                }
                
                #If the "From" of the edge is currently present, append to the previous entry.
                if (length(check) == 1)
                {
                        paths = c(paths, paste0(paths[check],":", edges[i,2]))
                }
        }
        last_node = unlist(lapply(strsplit(paths, ":"), function(x) x[length(x)]))
        paths = data.frame(Path = paths, Last = last_node)
        
        return(paths)
}

#This function takes the output of the parse.dot function, and the boolean table for the subclones, and returns a dataframe
#that may be used to annotate a seurat object as metadata with all the evolutionary branch information.
get.cells = function(parsed_dot, subclone_table)
{
        parsed_dot[,1] = gsub("RPL5G", "RPL5_G", parsed_dot[,1])
        parsed_dot[,1] = gsub("RPL5A", "RPL5_A", parsed_dot[,1])
        parsed_dot[,2] = gsub("RPL5G", "RPL5_G", parsed_dot[,2])
        parsed_dot[,2] = gsub("RPL5A", "RPL5_A", parsed_dot[,2])
        
        #For each last node entry, the corresponding cells for that node must be positive for ALL of the genes in the path.
        cell_memberships = list(length = nrow(parsed_dot))
        for (i in 1:nrow(parsed_dot))
        {
                last = as.character(parsed_dot[i,2])
                path = as.character(parsed_dot[i,1])
                if (last != "Unknown")
                {
                        path = unlist(strsplit(path, ":"))
                        #All cells by definition, are positive for "germline."
                        path = path[-1]
                        
                        cells_check = data.frame(subclone_table[,colnames(subclone_table) %in% path])
                        rownames(cells_check) = rownames(subclone_table)
                        cells = rownames(cells_check)[rowSums(cells_check) == ncol(cells_check)]
                        cell_memberships [[i]] = cells
                }
                
        }
        subclone_metadata = data.frame(matrix(nrow = nrow(subclone_table), ncol = length(cell_memberships)))
        rownames(subclone_metadata) = rownames(subclone_table)
        colnames(subclone_metadata) = parsed_dot$Path
        subclone_metadata[is.na(subclone_metadata)] = 0
        
        for (i in 1:length(cell_memberships))
        {
                index = rownames(subclone_metadata) %in% cell_memberships[[i]]
                subclone_metadata[index,i] = 1
        }
        #For checking purposes.
        #data.frame(Last = parsed.dot$Last, n = apply(subclone.metadata, 2, sum))
        return(subclone_metadata)
}

#Function for visualization of gene signature information projected onto subclone information (graphviz) in a phylogenetic tree.
#input of this function includes dot (the original .gv file), output of calc.signature.scores.list (signature_scores), and output of get.cells (subclone_metadata).

viz.scores = function(dot, signature_scores, subclone_metadata) 
{
        
        #first, compare the rownames of signature_scores and the rownames of subclone_metadata to make sure they line up.
        signature_scores = signature_scores[order(rownames(signature_scores)),]
        subclone_metadata = subclone_metadata[order(rownames(subclone_metadata)),]
        
        if (sum(rownames(subclone_metadata) == rownames(signature_scores)) != nrow(subclone_metadata))
        {
                print("Check the cell names included in the input to this function.")
                return(0)
        }
        
        #For each of the nodes, we want to get a mean signature of all the cells of that node.
        nodes = colnames(subclone_metadata)
        
        subclones_signatures_means = list()
        
        for (i in 1:length(nodes))
        {
                temp = signature_scores[as.logical(subclone_metadata [,i]),]
                subclones_signatures_means[[i]] = colMeans(temp)
        }
        subclones_signatures_means = t(as.data.frame(subclones_signatures_means))
        rownames(subclones_signatures_means) = colnames(subclone_metadata)
        
        #remove NA, convert to colour.
        subclones_signatures_means = subclones_signatures_means[!is.na(subclones_signatures_means[,1]),]
        subclones_signatures_means_color = apply(subclones_signatures_means, 2, function (x) number.to.colors(x,
                                                                                                              colors = c("white", "blue"), num = 100))
        rownames(subclones_signatures_means_color) = rownames(subclones_signatures_means)
        
        #Annotate in the colours into new .gv files.
        
        lines = rownames(subclones_signatures_means_color)
        lines = unlist(lapply(strsplit(lines, ":"), tail, n=1L))
        
        filenames = vector()
        for (k in 1:length(colnames(signature_scores)))
        { 
                new_dot = dot
                for (i in 1:length(lines))
                {
                        b = grep(lines[i], dot)
                        colour = subclones_signatures_means_color[i,k]
                        test = gsub("];", ", style = filled, fillcolor = ", dot[b])
                        test = paste0(test, "\"" ,colour, "\"",  "];")
                        new_dot[b] = test
                }
                filenames[k] =  paste0("new_protocol_full_mutations_only_", colnames(signature_scores)[k], ".gv")
                writeLines(new_dot, filenames[k])
        }
        return (filenames)
}
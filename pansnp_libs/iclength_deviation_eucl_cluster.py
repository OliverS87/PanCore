# Write the Rcscript required for intracluster length deviation clustering into the
# output folder
from os import path
from os import remove

class IclengthClusterRscript:
    def __init__(self, out_path):
        self.out_path = out_path

    def write_script(self):
        with open(path.join(self.out_path, "iclength_deviation_eucl_cluster.r"), "w") as out_f:
            for line in self.get_script():
                out_f.write(line+"\n")

    def remove_script(self):
        remove(path.join(self.out_path, "iclength_deviation_eucl_cluster.r"))

    # sed -e 's/.*/"&",/' iclength_deviation_eucl_cluster.r
    def get_script(self):
        script =[
            "#!/usr/bin/env Rscript",
            "suppressMessages(library(vegan))",
            "suppressMessages(library(reshape))",
            "suppressMessages(library(gplots))",
            "suppressMessages(library(RColorBrewer))",
            "suppressMessages(library(gtools))",
            "argv = commandArgs(trailingOnly=TRUE)",
            "if (length(argv)!=5) {",
            "  stop('Invalid nr. of cmdline args')",
            "} ",
            "",
            "ic_length_vector_path <- argv[1]",
            "out_file <- argv[2]",
            "out_hm <- argv[3]",
            "max_cluster_per_iteration <- as.numeric(argv[4])",
            "prefix <- argv[5]",
            "# Parse intracluster region length vector",
            "ic_len_data = data.frame(read.csv2(header=TRUE, sep = ',', file = ic_length_vector_path), stringsAsFactors = F)",
            "ic_len_data$IC.ID <- paste(ic_len_data$CB.L,ic_len_data$CB.R, sep = '.')",
            "ic_len_data.cast <- cast(ic_len_data, SI~IC.ID, value = 'LEN.IC', fun.aggregate = max)",
            "# Make the SI the row names",
            "rownames(ic_len_data.cast) <- ic_len_data.cast[,1]",
            "# and remove the first data column with the SIs",
            "# So that we don't cluster on the SI",
            "ic_len_data.cast[,1] <- NULL",
            "# We do not need to include the reference in the clustering",
            "ic_len_data.cast <- ic_len_data.cast[-1,]",
            "# Remove columns/IC regions where all lengths are the same -> No additional inforamtion",
            "ic_len_data.cast.filter <- Filter(function(x) length(unique(x))>1, ic_len_data.cast)",
            "# Convert data into a matrix",
            "ic_len_data.mtrx <- data.matrix(ic_len_data.cast.filter)",
            "rownames(ic_len_data.mtrx) <- rownames(ic_len_data.cast.filter)",
            "colnames(ic_len_data.mtrx) <- colnames(ic_len_data.cast.filter)",
            "",
            "# Order columns with mixedorder",
            "ic_len_data.mtrx <- ic_len_data.mtrx[,mixedorder(colnames(ic_len_data.mtrx))]",
            "",
            "# Plot heatmap",
            "if(out_hm!='NA'){",
            "my_palette <- colorRampPalette(c('black', 'lightgrey', 'cyan'))(n = 1000)",
            "png(filename = out_hm, width = 1024, height = 768)",
            "heatmap.2(ic_len_data.mtrx, distfun = function(x) vegdist(x, method='euclidian', na.rm = TRUE), ",
            "          hclustfun = function(x) hclust(x, method = 'ward.D2'), Colv = NA, key.title = 'Z-Plot', ",
            "          xlab = 'Intracore region', ylab = 'Sequence ID', ",
            "          main = paste(prefix, 'between core blocks length deviation', sep='\n'), na.rm = TRUE, key = T, ",
            "          col =my_palette, dendrogram = 'row',scale = 'column', trace = 'none')",
            "dev.off()}",
            "# Cluster: Euclidean distance, ward.D2 hierarchical clustering",
            "ic_len_data.mtrx.dist <-vegdist(ic_len_data.mtrx,method='euclidian', na.rm = TRUE)",
            "ic_len_data.mtrx.clust <- hclust(ic_len_data.mtrx.dist, method = 'ward.D2')",
            "",
            "# Cut the tree into cluster",
            "# The aim is to find clusters with more than one element, ",
            "# at max to max_cluster_per_iteration multiclusters",
            "# If that is not possible, go with the lowest k_val that creates ",
            "# the max. number of clusters with more than one element",
            "multiclstr_count <- c()",
            "for(k_val in c(2:nrow(ic_len_data.cast.filter))){",
            "  # Cut tree into k_val groups",
            "  print(k_val)",
            "  cluster <- cutree(ic_len_data.mtrx.clust, k=k_val)",
            "  # Count the number of groups with more than one element",
            "  count <- sum(table(cluster)>1)",
            "  # Add this count to the list",
            "  multiclstr_count[k_val] <- count",
            "  # If this count equals the max. number of groups we are looking for, stop here",
            "  if (count==max_cluster_per_iteration) break}",
            "# Get the k_val for which we achieved the maximum number of multi-element groups",
            "max_k <- which(multiclstr_count == max(multiclstr_count, na.rm = T))[1]",
            "# Cut the tree into max_k groups",
            "cluster <- cutree(ic_len_data.mtrx.clust, k=max_k)",
            "write.table(cluster, out_file, quote = FALSE)"
        ]
        return script

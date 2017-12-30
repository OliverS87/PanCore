# Write the Rcscript required for rearrangement clustering into the
# output folder
from os import path
from os import remove


class RearrangementJacCluster:
    def __init__(self, out_path):
        self.out_path = out_path

    def write_script(self):
        with open(path.join(self.out_path, "rearrangement_jac_cluster.r"), "w") as out_f:
            for line in self.get_script():
                out_f.write(line+"\n")

    def remove_script(self):
        remove(path.join(self.out_path, "rearrangement_jac_cluster.r"))

    # sed -e 's/.*/"&",/' rearrangement_jac_cluster.r
    def get_script(self):
        script = [
            "#!/usr/bin/env Rscript",
            "suppressMessages(library(vegan))",
            "suppressMessages(library(gtools))",
            "suppressMessages(library(gplots))",
            "suppressMessages(library(RColorBrewer))",
            "argv = commandArgs(trailingOnly=TRUE)",
            "if (length(argv)!=5) {",
            "  stop('Invalid nr. of cmdline args')",
            "} ",
            "",
            "feat_del_path <- argv[1]",
            "out_file <- argv[2]",
            "out_plot <- argv[3]",
            "max_cluster_per_iteration <- as.numeric(argv[4])",
            "prefix <- argv[5]",
            "",
            "# Parse deletion feature vector",
            "feat_del <- data.frame(read.csv2(header=TRUE, sep = ',', file = feat_del_path, row.names = 1))",
            "# We do not need to include the reference in the clustering",
            "feat_del <- feat_del[-1,]",
            "# Compute the jaccard distance",
            "feat_del.dist <-vegdist(feat_del,method='jaccard')",
            "# Hierarchical clustering",
            "feat_del.clust <- hclust(feat_del.dist, method = 'ward.D2')",
            "",
            "# Order columns with mixedorder",
            "feat_del.mtrx <- data.matrix(feat_del)",
            "feat_del.mtrx <- feat_del.mtrx[,mixedorder(colnames(feat_del.mtrx))]",
            "# Plot heatmap with deletions",
            "if(out_plot!='NA'){",
            "  png(filename = out_plot, width = 1024, height = 768)",
            "  my_palette <- colorRampPalette(c('black', 'lightgrey', 'lightgrey'))(n = 2)",
            "  plot <- heatmap.2(feat_del.mtrx, distfun = function(x) vegdist(x, method='jaccard', na.rm = TRUE), ",
            "                    hclustfun = function(x) hclust(x, method = 'ward.D2'), Colv = NA, key.title = 'Z-Plot', ",
            "                    xlab = 'Intracore region', ylab = 'Sequence ID', ",
            "                    main = paste(prefix, 'intra-core region deletions', sep='\n'), na.rm = TRUE, key = T, ",
            "                    col =my_palette, dendrogram = 'row',scale = 'column', trace = 'none')",
            "  ",
            "  dev.off()}",
            "",
            "# Cut the tree into cluster",
            "# The aim is to find clusters with more than one element, at max to max_cluster_per_iteration multiclusters",
            "# If that is not possible, go with the lowest k_val that creates the max. number of clusters with more than one element",
            "multiclstr_count <- c()",
            "for(k_val in c(2:nrow(feat_del))){",
            "  # Cut tree into k_val groups",
            "  cluster <- cutree(feat_del.clust, k=k_val)",
            "  # Count the number of groups with more than one element",
            "  count <- sum(table(cluster)>1)",
            "  # Add this count to the list",
            "  multiclstr_count[k_val] <- count",
            "  # If this count equals the max. number of groups we are looking for, stop here",
            "  if (count==max_cluster_per_iteration) break}",
            "# Get the k_val for which we achieved the maximum number of multi-element groups",
            "max_k <- which(multiclstr_count == max(multiclstr_count, na.rm = T))[1]",
            "# Cut the tree into max_k groups",
            "cluster <- cutree(feat_del.clust, k=max_k)",
            "write.table(cluster, out_file, quote = FALSE)"
        ]
        return script

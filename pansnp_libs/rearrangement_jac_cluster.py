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
            "argv = commandArgs(trailingOnly=TRUE)",
            "if (length(argv)!=5) {",
            "  stop('Invalid nr. of cmdline args')",
            "} ",
            "",
            "feat_del_path <- argv[1]",
            "out_file <- argv[2]",
            "out_plot <- argv[3]",
            "min_cluster_goal <- as.numeric(argv[4])",
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
            "# Plot tree",
            "if(out_plot!='NA'){",
            "  png(filename = out_plot, width = 1024, height = 768)",
            "  plot(as.dendrogram(feat_del.clust), xlab = 'Sequence ID (SI)', ylab = 'Jaccard distance', main = paste('Rearrangement clustering for\n', prefix, '\nusing', feat_del.clust$method))",
            "  dev.off()}",
            "# Cut the tree into cluster",
            "# Try to get at min. two clusters with more than one member",
            "for(k_val in c(1:nrow(feat_del))){",
            "  cluster <- cutree(feat_del.clust, k=k_val)",
            "  if (sum(table(cluster)>1)>=min_cluster_goal) break",
            "}"
        ]
        return script

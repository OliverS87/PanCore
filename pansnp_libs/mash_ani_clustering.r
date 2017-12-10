#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly=TRUE)
if (length(argv)!=5) {
  stop('Invalid nr. of cmdline args')
} 

mash_dist_path <- argv[1]
out_file <- argv[2]
out_tree <- argv[3]
max_cluster_per_iteration <- as.numeric(argv[4])
prefix <- argv[5]

# mash to tree
mash_data = as.data.frame(read.csv2(mash_dist_path, sep='', stringsAsFactors = F))
# Make the SI the row names
rownames(mash_data) <- mash_data[,1]
# and remove the first data column with the SIs
# So that we don't cluster on the SI
mash_data[,1] <- NULL
colnames(mash_data) <- rownames(mash_data)
mash_data.m <- as.matrix(mash_data)
mash_data.d <- as.dist(mash_data.m)
# Cluster: Euclidean distance, ward.D2 hierarchical clustering
mash_data.c <- hclust(mash_data.d, method = 'ward.D2')
# Plot heatmap
png(filename = out_tree, width = 1024, height = 768)
plot(as.dendrogram(mash_data.c), xlab = 'Sequence ID (SI)', ylab = 'Distance', main = paste('ANI clustering for\n', prefix, '\nusing', mash_data.c$method))
dev.off()

# Cut the tree into cluster
# The aim is to find clusters with more than one element, at max to max_cluster_per_iteration multiclusters
# If that is not possible, go with the lowest k_val that creates the max. number of clusters with more than one element
multiclstr_count <- c()
for(k_val in c(2:nrow(mash_data.m))){
  # Cut tree into k_val groups
  cluster <- cutree(mash_data.c, k=k_val)
  # Count the number of groups with more than one element
  count <- sum(table(cluster)>1)
  # Add this count to the list
  multiclstr_count[k_val] <- count
  # If this count equals the max. number of groups we are looking for, stop here
  if (count==max_cluster_per_iteration) break}
# Get the k_val for which we achieved the maximum number of multi-element groups
max_k <- which(multiclstr_count == max(multiclstr_count, na.rm = T))[1]
# Cut the tree into max_k groups
cluster <- cutree(mash_data.c, k=max_k)
write.table(cluster, out_file, quote = FALSE)

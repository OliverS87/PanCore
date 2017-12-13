# Randomize a set of clusters
# Random clusters are used to assess the significance
# of computed clusters
from random import shuffle
class RandomCluster:
    def randomize(self, cluster_list):
        shuffled_files = [file for cluster in cluster_list for file in cluster]
        shuffle(shuffled_files)
        cluster_lenghts = [len(lst) for lst in cluster_list]
        shuffled_cluster = []
        next_pos = 0
        for length in cluster_lenghts:
            shuffled_cluster.append(shuffled_files[next_pos:next_pos+length])
            next_pos += length
        return shuffled_cluster

if __name__ == '__main__':
    print("welcome to maine")
    rc = RandomCluster()
    my_clstr = [["a"], ["d", "e"], ["f"]]
    my_r_clstr = rc.randomize(my_clstr)
    print(my_clstr)
    print(my_r_clstr)




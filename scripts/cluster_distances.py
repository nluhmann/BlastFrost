from ete2 import Tree
import sys
import collections
import itertools

# read tree file
file = open(sys.argv[1], "r")
newick = file.read()
file.close()
tree = Tree(newick, format=1)
#print tree.get_ascii()
#print " "


for node in tree.traverse(strategy="postorder"):
    if node.is_leaf():
        node.add_feature("clusterID",0)




# read clustering file, then annotate each leaf in the tree by its clusterID

clusters = collections.defaultdict(list)
f = open(sys.argv[2],"r")
for line in f:
    arr = line.rstrip("\n").split("\t")
    id_arr = arr[0].split("/")[-1].split("_")
    #print id_arr
    strain_id = id_arr[0]+"_"+id_arr[1]+"_AS"
    #print strain_id
    cluster_id = arr[2]
    leaf = tree.get_leaves_by_name(strain_id)
    if len(leaf) > 1 :
        print "ERROR"
    leaf[0].add_feature("clusterID",cluster_id)
    clusters[cluster_id].append(strain_id)
    #print leaf[0].clusterID

f.close()



#for cluster in clusters:
#    if (len(clusters[cluster]) == 1):
#        print clusters[cluster]
#        print tree.check_monophyly(clusters[cluster], target_attr="name")
#        if (tree.check_monophyly(clusters[cluster], target_attr="name")[0]):
#            print cluster



#compute max. intra-cluster distance

f = open("intra_cluster_dist","w")
for cluster in clusters:
    max_dist = 0
    for pair in itertools.combinations(clusters[cluster],2):
        a = tree&pair[0]
        dist = a.get_distance(pair[1])
        if (dist > max_dist):
            max_dist = dist

    f.write(str(cluster)+"\t"+str(max_dist)+"\t"+str(max_dist/len(clusters[cluster]))+"\n")
f.close()
            

f = open(sys.argv[3],"r")
distance = {}
for line in f:
            arr = line.rstrip("\n").split("\t")
            pair = tuple(arr[0:2])
            distance[pair] = arr[2]
f.close()


            
f = open("inter_cluster_dist","w")
for pair in itertools.combinations(clusters.keys(),2):
    #find min distance between these two
    max_dist = 0
    if pair in distance:
        hamming = distance[pair]
    elif pair[::-1] in distance:
        hamming = distance[(pair[1],pair[0])]
    else:
        print "ERROR"
    for it in itertools.product(clusters[pair[0]], clusters[pair[1]]):
        a = tree&it[0]
        dist = a.get_distance(it[1])
        if (max_dist == 0 or dist < max_dist):
            max_dist = dist

    f.write(pair[0]+"\t"+pair[1]+"\t"+str(max_dist)+"\t"+str(hamming)+"\t"+str(max_dist/int(hamming))+"\n")
f.close()

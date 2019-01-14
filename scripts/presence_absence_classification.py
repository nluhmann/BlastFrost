import sys
import collections
import itertools

#for each input file, save list of present strains


markers = collections.defaultdict(list)
for file in sys.argv[1:]:
    f = open(file, "r")
    name = file.split("_")[-1]
    for line in f:
        strain = line.split("\t")[1]
        markers[name].append(strain)


    f.close()


all_strains = {}

for elem in markers:
    for strain in markers[elem]:
        if not strain in all_strains:
            all_strains[strain] = []


order = []

for elem in markers:
    order.append(elem)
    for strain in all_strains:
        if strain in markers[elem]:
            all_strains[strain].append('1')
        else:
            all_strains[strain].append('0')




out = open("outfile.txt","w")
out.write("strain" + "\t"+"\t".join(order)+"\n")
for elem in all_strains.keys():
    out.write(elem+"\t"+"\t".join(all_strains[elem])+"\n")

out.close()
            


clusters = collections.defaultdict(list)

for elem in all_strains:
    clusters[tuple(all_strains[elem])].append(elem)

clusterIDs = {}

out = open("clusters.txt","w")
clusterID = 0
for elem in clusters:
    for strain in clusters[elem]:
        out.write(strain + "\t" + str(elem) + "\t" + str(clusterID)+"\n")
    clusterIDs[elem] = clusterID
    clusterID+=1

out.close()


out = open("cluster_dist.txt","w")

#ToDo: I need the cluster ID for this, not the cluster element....
for pair in itertools.combinations(clusters.keys(),2):
    dist = 0
    for i in range(0,len(pair[0])):
        if pair[0][i] != pair[1][i]:
            dist += 1

    out.write(str(clusterIDs[pair[0]]) + "\t" + str(clusterIDs[pair[1]])+"\t"+str(dist)+"\n")

out.close()






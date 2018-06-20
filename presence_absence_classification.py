import sys
import collections


#for each input file, save list of present strains


markers = collections.defaultdict(list)
for file in sys.argv[1:]:
    f = open(file, "r")
    name = file.split("_")[-1]

    for line in f:
        strain = line.split("\t")[0]
        markers[name].append(strain)


    f.close()


all_strains = {}

for elem in markers:
    for strain in markers[elem]:
        if not strain in all_strains:
            all_strains[strain] = []



for elem in markers:
    for strain in all_strains:
        if strain in markers[elem]:
            all_strains[strain].append('1')
        else:
            all_strains[strain].append('0')




out = open("outfile.txt","w")

for elem in all_strains.keys():
    out.write(elem+"\t"+"\t".join(all_strains[elem])+"\n")

out.close()
            


clusters = collections.defaultdict(list)

for elem in all_strains:
    clusters[tuple(all_strains[elem])].append(elem)


out = open("clusters.txt","w")
clusterID = 0
for elem in clusters:
    for strain in clusters[elem]:
        out.write(strain + "\t" + str(elem) + "\t" + str(clusterID)+"\n")

    clusterID+=1




out.close()



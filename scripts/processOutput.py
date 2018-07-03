import sys
import collections
import itertools

#read order of assembly barcodes in newick tree, assuming that subtrees are clustered together in the file format
with open(sys.argv[1], 'r') as nwkfile:
    tree = nwkfile.read().replace('\n', '')

leaf_order = []
arr = tree.split(",")

for elem in arr:
    x = elem.split(":")[0]
    if "(" in x:
        leaf = x.split("(")[-1]
    else:
        leaf = x
    leaf_order.append(leaf)

#read metadata table (to save: assembly barcode - 42, subspecies - 21, serovar - 75)
subspecies = {}
serovars = {}
serovar_clusters = collections.defaultdict(list)
f = open(sys.argv[3],"r")
for line in f:
    if not "Uberstrain" in line:
        arr = line.rstrip("\n").split("\t")
        if arr[20] == "":
            arr[20] = "NA"
        subspecies[arr[41]] = arr[20]
        if arr[74] == "":
            arr[74] = "NA"
        if "|" in arr[74]:
            s = arr[74].split("|")
            arr[74] = "|".join(s[0:2])
            print arr[74]
        serovars[arr[41]] = arr[74]
        if arr[74] != "NA":
            serovar_clusters[arr[74]].append(arr[41])
f.close()

#read gene names (extracted from header file)
geneNames = {}
queries = set()
f = open(sys.argv[4],"r")
for line in f:
    arr = line.split(" ")
    gi = arr[0][1:]
    gene = arr[1]
    geneNames[gi]=gene
    queries.add(gi)
f.close()

#read search results
serovar_constistency = {}
search_result = {}
f = open(sys.argv[2],'r')
for line in f:
    arr = line.rstrip("\n").split("\t")
    strain_arr = arr[1].split("/")[-1].split("_")
    strain = "_".join(strain_arr[0:3])
    subsp = subspecies[strain]
    serov = serovars[strain]
    gi = arr[0]
    gene = geneNames[gi]
    line = gi+"\t"+strain+"\t"+arr[2]+"\t"+arr[3]+"\t"+subsp+"\t"+gene
    search_result[(strain,gi)] = line
f.close()

#extend search_result to include everything not found
#serovar_results = collections.defaultdict(list)
#for leaf in leaf_order:
#    for gi in queries:
#        if (leaf,gi) in search_result:
#            present = 1
#            line = search_result[(leaf,gi)]
#            line+="\t1\n"
#        else:
#            present = 0
#            line = gi+"\t"+leaf+"\t-\t-\t"+subspecies[leaf]+"\t"+geneNames[gi]+"\t0\n"
#        serovar_results[serovars[leaf]].append(line)


#translate leaf_order into serovar order




#write serovar clusters in order, keeping one representative per cluster


out = open(sys.argv[2]+"_processed","w")
#out.write("query\tstrain\tevalue\tkmer_hits\tsubspecies\tgene\tpresence\n")
query_labels = []
for gi in queries:
    gene = geneNames[gi]
    new = gi.split("(")[1][:-1]+" "+gene
    query_labels.append(new)
out.write("strain\t"+"\t".join(query_labels)+"\tsubspecies\tserovar\n")
for leaf in leaf_order:
    subspec = subspecies[leaf]
    sov = serovars[leaf]
    out.write(leaf+"\t")
    for gi in queries:
        if (leaf,gi) in search_result:
            #line = search_result[(leaf,gi)]
            out.write("1\t")
        else:
            #no hit
            #out.write(gi+"\t"+leaf+"\t-\t-\t"+subspecies[leaf]+"\t"+geneNames[gi]+"\t0\n")
            out.write("0\t")
    out.write(subspec+"\t"+sov+"\n")
out.close()

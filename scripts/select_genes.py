import sys



#read gene header to extract
header = []

f = open(sys.argv[1],"r")
for line in f:
    header.append(line.rstrip("\n"))

f.close()


hits = {}
seq = ""
f = open(sys.argv[2],"r")
for line in f:
    if line.startswith(">"):
        if seq != "":
            if head in header:
                hits[head] = seq
            seq = ""
        head = line.rstrip("\n")
    else:
        seq += line.rstrip("\n")


f.close()


out = open(sys.argv[3],"w")
for hit in hits:
    out.write(hit+"\n")
    out.write(hits[hit]+"\n")

out.close()
    

import sys
import random
import collections

#parameters
k = 13
sigma = 4
len_x = 20000
len_y = 20000




def random_kmers(a,k):
    nucl = ['A','C','G','T']

    kmer_hash = collections.defaultdict(int)

    for x in range(0, a):
        #print x
        new_kmer = ""
        for y in range(0,k):
            #draw random nucleotide
            new_kmer += random.choice(nucl)

        kmer_hash[new_kmer] += 1
        #print new_kmer

    return kmer_hash



kmer_x = random_kmers(len_x,k)
kmer_y = random_kmers(len_y,k)

matches = 0


keys_x = kmer_x.keys()

for kmer in kmer_y.keys():
    if kmer in keys_x:
        matches+=1



print "k: "+str(k)+"\n"
print "len_x: "+str(len_x)+"\n"
print "len_y: "+str(len_y)+"\n"
print "matches: "+str(matches)+"\n"











# BlastFrost

BlastFrost is a highly efficient method for querying 100,000s of genome assemblies. BlastFrost builds on the recently developed Bifrost, which generates a dynamic data structure for compacted and colored de Bruijn graphs from bacterial genomes. BlastFrost queries a Bifrost data structure for sequences of interest, and extracts local subgraphs, thereby enabling the efficient identification of the presence or absence of individual genes or single nucleotide sequence variants.

You can learn more about Bifrost graphs and how to build them [here](https://github.com/pmelsted/bifrost)

## Requirements
* precomputed Bifrost graph as input
* query sequences of interest
* C++11 compiler


## Install from source
1. Clone the github repository
```
git clone --recursive https://github.com/nluhmann/BlastFrost.git
cd BlastFrost
```

2. Install BlastFrost
```
mkdir build; cd build; cmake ..; make
```

Note that this will also compile Bifrost(https://github.com/pmelsted/bifrost) for you. You can use the executable in BlastFrost/bifrost/src/Bifrost to create your bifrost graph that is required as input to Bifrost.

## Conda
TBA


### Large k-mers
Bifrost has a standard maximum k-mer size of 31, but can be adjusted at compilation to accept larger values for k. If you compiled Bifrost with a larger max k, you will accordingly need to adjust the parameter in CMakeList.txt:

```
set(-MAX_KMER_SIZE=64)
```
In this example, you can then use BlastFrost for graphs with a maximum k-mer size of 63. 

### Install Troubleshooting
You might run into problems linking BlastFrost to your Bifrost installation (you are probably running MacOS?).
Please see the [Bifrost troubleshooting section](https://github.com/pmelsted/bifrost#troubleshooting) to correctly set environment variables on your system.



## Usage
```
BlastFrost -g <BifrostGraph> -f <BifrostColors> -q <query_sequences> -o <outfile_prefix>

Mandatory parameters with required argument:

  -g,         Input Bifrost graph file (GFA format)
  -f,         Input Bifrost color file (BFG_COLORS format)
  -q,         Query sequences (multiple FASTA)
  -o,         Prefix for search result files (default: 'output')

Optional parameters:

  -e,          Enable subgraph extraction
  -t,          Number of threads (default is 1)
  -d,          Compute and search d-neighborhood of queried k-mers
  -c,          Enhance gfa file with color information, option can be run without a query search.
  -s,          Average size of genomes in Bifrost graph in Mb
  -a,          Return file containing all colors currently present in Bifrost graph
  -v,          Print information messages during construction

```



## Example commands
For testing, the directory `testdata` in this repository contains a pre-computed Bifrost graph with *k = 31* and a query sequence in fasta format. 

### K-mer search
Simple k-mer search for query on a single thread, providing the average genome length of *Y. pestis* as input:
```
/path/to/BlastFrost -g yersinia_5.gfa -f yersinia_5.bfg_colors -q testquery -o outtest -s 4.6 -v -t 1
```

### Subgraph extraction
We can extend the k-mer search result by extracting a subgraph for the query, additionally writing path sequences to ouput:
```
/path/to/BlastFrost -g testdata/yersinia_cdBG.gfa -f testdata/yersinia_cdBG.bfg_colors -q testdata/testquery -o outtest -t 1 -s 4.6 -v -e
```

## BlastFrost-AMR
BlastFrost-AMR is an application of BlastFrost, based on AMR genes in [Bacterial Antimicrobial Resistance Reference Gene Database](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047). 

BlastFrost-AMR can be run directly in any Python 3 environment. First, download the lastest database from NCBI:  
```
BlastFrost-AMR build -b </path/to/bifrost>
```
 

To search the presence/absence of AMR genes in a pre-computed Bifrost graph: 
```
BlastFrost-AMR query -g </path/to/graph> -p <prefix>
```

Also, use '--help' to get detailed manual for BlastFrost-AMR. 


## Citation
[![DOI](https://zenodo.org/badge/135275732.svg)](https://zenodo.org/badge/latestdoi/135275732)

Our BlastFrost paper is currently under review, but you can check out the preprint:

[BlastFrost: Fast querying of 100,000 s of bacterial genomes in Bifrost graphs](https://www.biorxiv.org/content/10.1101/2020.01.21.914168v1.abstract)

[Bifrost - Highly parallel construction and indexing of colored and compacted de Bruijn graphs](https://www.biorxiv.org/content/10.1101/695338v2.abstract)

```
@article{luhmann2020blastfrost,
  title={BlastFrost: Fast querying of 100,000 s of bacterial genomes in Bifrost graphs},
  author={Luhmann, Nina and Holley, Guillaume and Achtman, Mark},
  journal={BioRxiv},
  year={2020},
  publisher={Cold Spring Harbor Laboratory}
}

@article{holley2020bifrost,
  title={Bifrost: highly parallel construction and indexing of colored and compacted de Bruijn graphs},
  author={Holley, Guillaume and Melsted, P{\'a}ll},
  journal={Genome biology},
  volume={21},
  number={1},
  pages={1--20},
  year={2020},
  publisher={BioMed Central}
}
```


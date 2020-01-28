# BlastFrost

BlastFrost is a highly efficient method for querying 100,000s of genome assemblies. BlastFrost builds on the recently developed Bifrost, which generates a dynamic data structure for compacted and colored de Bruijn graphs from bacterial genomes. BlastFrost queries a Bifrost data structure for sequences of interest, and extracts local subgraphs, thereby enabling the efficient identification of the presence or absence of individual genes or single nucleotide sequence variants.

You can learn more about Bifrost graphs and how to build them [here](https://github.com/pmelsted/bifrost)

## Requirements
* precomputed Bifrost graph as input
* query sequences of interest
* C++11 compiler


## Install
First, you need to install Bifrost: https://github.com/pmelsted/bifrost

Then you can install BlastFrost:

```
mkdir build; cd build; cmake ..; make
```



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
  -k,          Length of k-mers (default is 31, max. is 31)
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
/path/to/BlastFrost -g testdata/yersinia_cdBG.gfa -f testdata/yersinia_cdBG.bfg_colors -q testdata/testquery -o outtest -t 1 -s 4.6 -v 
```

### Subgraph extraction
We can extend the k-mer search result by extracting a subgraph for the query, additionally writing path sequences to ouput:
```
/path/to/BlastFrost -g testdata/yersinia_cdBG.gfa -f testdata/yersinia_cdBG.bfg_colors -q testdata/testquery -o outtest -t 1 -s 4.6 -v -e
```





## Citation




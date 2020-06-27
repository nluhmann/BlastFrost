#!/usr/bin/env python
import sys, os, logging, subprocess
from pip._internal import main as pipmain
try :
    import click, requests
except :
    pipmain(['install', 'click', 'request'])
    import click, requests

logging.basicConfig(format='%(asctime)s | %(message)s', level='INFO')

dirname = os.path.dirname(__file__)
blastFrost = os.path.join(dirname, 'build', 'BlastFrost')
barrd = os.path.join(dirname, 'barrd')

def checkFiles() :
    if not os.path.isfile(blastFrost):
        logging.error('BlastFrost has not been compiled. Please compile it before running BlastFrost-AMR')
        sys.exit(1)
    if not os.path.isfile(barrd+'.nuc') or not os.path.isfile(barrd+'.incompatible') :
        logging.error('Bacterial Antimicrobial Resistance Reference Gene Database is not found. Use BlastFrost-AMR -b <path/to/bifrost> to build the database.')
        sys.exit(1)

@click.group()
def main() :
    pass

@main.command()	
@click.option('-b', '--bifrost', help='path to Bifrost', default='Bifrost')
@click.option('-t', '--thread', help='[Default: 1] Number of threads. defulat: 1', default=1, type=int)
def build(bifrost, thread) :
    logging.warning('Downloading database from from https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS')
    with open(barrd+'.nuc', 'wt') as fout :
        fout.write(requests.get('https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS').text)
    logging.warning('Building database')
    if not os.path.isdir(barrd+'.dir') :
        os.makedirs(barrd+'.dir')
    fout, names = None, []
    with open(barrd+'.nuc', 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split('|')[0]
                sname = os.path.join(barrd+'.dir', name)
                names.append( sname )
                if fout != None :
                    fout.close()
                fout = open(sname, 'wt')
            fout.write(line)
    fout.close()
    with open(barrd+'.txt', 'wt') as fout :
        fout.write('\n'.join(names))

    subprocess.Popen('{0} build -r {1} -k 21 -c -o {2} -t {3}'.format(
        bifrost, barrd+'.txt', barrd, thread
    ).split(), stdout=subprocess.PIPE).communicate()
    subprocess.Popen('{0} -g {1} -f {2} -q {3} -d 0 -o {4} -t {5}'.format(
        blastFrost, barrd+'.gfa', barrd+'.bfg_colors', barrd+'.nuc', barrd, thread
    ).split(), ).communicate()
    conflict_pair = {}
    with open(barrd+'_barrd.nuc.search', 'rt') as fin :
        for line in fin :
            part = line.strip().split('\t')
            q, r = part[0].split('|')[0], part[1].rsplit('/', 1)[-1]
            if q ==r :
                continue
            if q not in  conflict_pair :
                conflict_pair[q] = {}
            if r not in conflict_pair :
                conflict_pair[r] = {}
            conflict_pair[q][r] = 1
            conflict_pair[r][q] = 1
    
    with open(barrd+'.incompatible', 'wt') as fout :
        for k, y in conflict_pair.items() :
            fout.write('{0}\t{1}\n'.format(k, ','.join(y)))
    for fn in (barrd+'.txt', barrd+'.gfa', barrd+'.bfg_colors', barrd+'_barrd.nuc.search') :
        os.unlink(fn)
    import shutil
    shutil.rmtree(barrd+'.dir')
    
    
@main.command()
@click.option('-g', '--graph', help='[Required] Genomic graph built with "bifrost -c"')
@click.option('-p', '--prefix', help='[Optional; Default: same as graph] prefix of the output', default='')
@click.option('-d', '--dist', help='[Optional; Default: 1] "-d" parameter for BlastFrost. default: 1', default=1, type=int)
@click.option('-n', '--nofilter', help='[Default]: filter out similar hits; Flag: no filter', is_flag=True, default=False)
@click.option('-t', '--thread', help='[Default: 1] Number of threads. defulat: 1', default=1, type=int)
def query(graph, prefix, dist, nofilter, thread) :
    checkFiles()
    if prefix == '' :        
        prefix = graph
        
    incompatible = {}
    with open(barrd+'.incompatible', 'rt') as fin :
        for line in fin :
            k, v = line.strip().split('\t')
            incompatible[k] = v.split(',')

    colors = graph[:-3] + 'bfg_colors'
    logging.info('BlastFrost starts...')
    subprocess.Popen('{0} -g {1} -f {2} -q {3} -o {4} -d {5} -t {6}'.format(
        blastFrost, graph, colors, barrd+'.nuc', prefix, dist, thread
        ).split()).communicate()

    logging.info('BlastFrost finishes. Parsing the results...')
    results = []
    with open('{0}_{1}.nuc.search'.format(prefix, os.path.basename(barrd)), 'r') as fin :
        for line in fin :
            gene, genome, evalue, aln = line.strip().split('\t')
            score = 0
            alns = [[int(t) for t in a.split(':')] for a in aln[:-1].split(',')]
            score = sum([(3-t)*l for t, l in alns if t > 0])
            size = sum([ l for t, l in alns ])
            if score < size :
                continue

            gene_info = gene.split('|')
            allele, family, desc = gene_info[5:8]
            results.append([genome, -score, family, allele, gene_info[0], evalue, desc, aln])
    nr = {}
    results.sort()
    logging.info('Writing outputs into {0}...'.format(prefix+'.BARRD.result'))    
    with open(prefix+'.BARRD.result', 'wt') as fout :
        fout.write('#Genome\tFamily\tAllele\tReference\tE-value\tMatches\tDescription\n')
        for genome, score, family, allele, ref, evalue, desc, aln in results :
            if (genome, ref) in nr or (genome, family) in nr :
                continue
            if not nofilter :
                nr[(genome, family)] = 1
                for k in incompatible.get(ref, {}) :
                    nr[(genome, k)] = 1

            fout.write('\t'.join([genome, family, allele, ref, evalue, aln, desc])+'\n')


if __name__ == '__main__' :
    main()

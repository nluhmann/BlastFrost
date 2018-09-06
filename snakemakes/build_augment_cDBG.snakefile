
#read all assembly files into list
#NOTE: the graph will be build on all these input assemblies
ASSEMBLIES, = glob_wildcards("assemblies/{id}.fasta")

#name the current run to diff output files
RUN_ID = "assembly1"

#set path for BlastFrost
BLASTFROST_PATH = "/home/nina/my_software/BlastFrost/build/BlastFrost"








rule all:
     input:
         expand("{RUN_ID}_joint_cDBG_colored.gfa", RUN_ID=RUN_ID)

	

rule build_cDBG:
     input:
        expand("assemblies/{assembly}.fasta",assembly=ASSEMBLIES)
     output:
        "{RUN_ID}_joint_cDBG.gfa"
     threads: 30
     benchmark:
         "benchmarks/Bifrost/{RUN_ID}_joint_cDBG.tsv"
     log:
         "logs/Bifrost/{RUN_ID}_joint_cDBG.log"
     message:
          "Run Bifrost with {threads} threads on all assemblies. See logfile {log}."
     shell:
         "Bifrost build -v -k 31 -c -t {threads} -o {RUN_ID}_joint_cDBG -r {input} > {log}"



rule augmentGFA:
    input:
        "{RUN_ID}_joint_cDBG.gfa"
    output:
        "{RUN_ID}_joint_cDBG_colored.gfa"
    threads: 30
    benchmark:
        "benchmarks/BlastFrost/{RUN_ID}_joint_cDBG_colored.tsv"
    log:
        "logs/BlastFrost/{RUN_ID}_joint_cDBG_colored.log"
    message:
        "Run BlastFrost with {threads} threads to augment gfa file with colors."
    shell:
        "{BLASTFROST_PATH} -v -t {threads} -g {RUN_ID}_joint_cDBG.gfa -f {RUN_ID}_joint_cDBG.bfg_colors -c"

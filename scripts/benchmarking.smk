"""
Dereplication benchmarking pipeline

pixi run snakemake \
    --snakefile scripts/benchmarking.smk \
    --directory results \
    --profile aqua \
    --keep-going --rerun-triggers mtime --cores 64 --local-cores 1
"""

import os
import polars as pl

GENOME_COUNTS = [100, 500, 1000, 5000, 10000, 50000]
REPEATS = list(range(1, 6))
TMPDIR = "/scratch/microbiome/aroneys/tmp"

# Genomes are sourced from https://doi.org/10.1038/s41592-025-02901-1
genomes = (
    pl.read_csv("/work/microbiome/ibis/SRA_genomes/compiled/genome_metadata.tsv", separator="\t")
    .filter(pl.col("completeness") - 5 * pl.col("contamination") >= 50)
    .select(
        genome = "name",
        fasta = pl.concat_str(pl.lit("/work/microbiome/ibis/SRA_genomes/"), pl.col("genome_path")),
        hash = pl.struct("name", "genome_path").hash(),
        )
    .sort("hash")
)

# GlobDB r232 reference genomes
GLOBDB_DIR = "/work/microbiome/db/globdb/232/globdb_r232_genome_fasta"
GLOBDB_CHUNK_DIRS = sorted(
    os.path.join(GLOBDB_DIR, d)
    for d in os.listdir(GLOBDB_DIR)
    if os.path.isdir(os.path.join(GLOBDB_DIR, d))
)

#################
### Functions ###
#################
def get_mem_mb(wildcards, threads, attempt):
    mem = 8 * 1000 * threads * attempt
    if mem == 1024000:
        return 1000000
    elif mem == 512000:
        return 500000
    elif mem == 256000:
        return 250000
    else:
        return mem

def get_genomes(wildcards):
    return (
        genomes
        .sample(fraction=1.0, shuffle=True, seed=int(wildcards.repeat))
        .head(int(wildcards.genome_count))
        .get_column("fasta")
        .to_list()
    )

####################
### Global rules ###
####################
# Compare results to skani with and without --slow
rule all:
    input:
        "galah/done/all_done",
        "galah_low_memory/done/all_done",
        "skder_greedy/done/all_done",
        "pyani_ANIm/done/all_done",
        "drep_skani/done/all_done",
        # ----------------------------------------------------------------------
        # Reference-based benchmarks: sample genomes clustered together with the
        # full GlobDB r232 reference set (~346k genomes).
        "galah_globdb/done/all_done",
        "galah_low_memory_globdb/done/all_done",
        "skder_greedy_globdb/done/all_done",
        # pyani_ANIm is skipped
        "drep_skani_globdb/done/all_done",
        "galah_reference/done/all_done",

rule compile_done:
    input:
        expand("{{tool}}/cluster_{genome_count}_rep{repeat}", genome_count=GENOME_COUNTS, repeat=REPEATS),
    localrule: True
    output:
        touch("{tool}/done/all_done"),

rule list_genomes:
    output:
        "genome_lists/genomes_{genome_count}_rep{repeat}.txt",
    localrule: True
    params:
        genomes = get_genomes,
    run:
        with open(output[0], "w") as f:
            for genome in params.genomes:
                f.write(f"{genome}\n")

rule link_genomes:
    output:
        directory("genome_links/genomes_{genome_count}_rep{repeat}"),
    localrule: True
    params:
        genomes = get_genomes,
    run:
        os.makedirs(output[0], exist_ok=True)
        for genome in params.genomes:
            genome_name = os.path.basename(genome)
            genome_link = os.path.join(output[0], genome_name)
            if not os.path.exists(genome_link):
                os.symlink(genome, genome_link)

rule list_globdb_genomes:
    output:
        "genome_lists/globdb_reference.txt",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 6*60*attempt,
    shell:
        "find {GLOBDB_DIR} -name '*.fa.gz' | sort > {output} "

rule list_genomes_with_reference:
    input:
        sample = "genome_lists/genomes_{genome_count}_rep{repeat}.txt",
        reference = "genome_lists/globdb_reference.txt",
    output:
        "genome_lists_with_reference/genomes_{genome_count}_rep{repeat}.txt",
    localrule: True
    shell:
        "cat {input.sample} {input.reference} > {output} "

###########################
### Dereplication tools ###
###########################
rule galah:
    input:
        "genome_lists/genomes_{genome_count}_rep{repeat}.txt",
    output:
        directory("galah/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/galah/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/galah/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e galah "
        "galah cluster "
        "--genome-fasta-list {input} "
        "--ani 95 "
        "--min-aligned-fraction 15 "
        "--skip-sanitise-headers "
        "--output-cluster-definition {output}/clusters.tsv "
        "--threads {threads} "
        "&> {log} "

rule galah_low_memory:
    input:
        "genome_lists/genomes_{genome_count}_rep{repeat}.txt",
    output:
        directory("galah_low_memory/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/galah_low_memory/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/galah_low_memory/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e galah "
        "galah cluster "
        "--genome-fasta-list {input} "
        "--ani 95 "
        "--min-aligned-fraction 15 "
        "--skip-sanitise-headers "
        "--output-cluster-definition {output}/clusters.tsv "
        "--threads {threads} "
        "--low-memory "
        "&> {log} "

rule skder_greedy:
    input:
        "genome_links/genomes_{genome_count}_rep{repeat}",
    output:
        directory("skder_greedy/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/skder_greedy/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/skder_greedy/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e skder "
        "skder "
        "--genomes {input} "
        "--percent-identity-cutoff 95 "
        "--aligned-fraction-cutoff 15 "
        "--dereplication-mode greedy "
        "--symlink "
        "--output-directory {output}/skder "
        "--threads {threads} "
        "&> {log} "

rule pyani_ANIm:
    input:
        "genome_links/genomes_{genome_count}_rep{repeat}",
    output:
        directory("pyani_ANIm/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/pyani_ANIm/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/pyani_ANIm/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e pyani "
        "pyani-plus anim "
        "$(pwd)/{input} "
        "--create-db "
        "--database {output}/pyani.db "
        # No threads arg? All by default?
        # "--threads {threads} "
        "&> {log} "
        "&& "
        "TMPDIR={TMPDIR} "
        "pixi run -e pyani "
        "pyani-plus classify "
        "--database {output}/pyani.db "
        "--outdir {output}/classify "
        "&>> {log} "

rule drep_skani:
    input:
        "genome_lists/genomes_{genome_count}_rep{repeat}.txt",
    output:
        directory("drep_skani/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/drep_skani/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/drep_skani/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e drep "
        "dRep dereplicate "
        "{output} "
        "--genomes {input} "
        "--ignoreGenomeQuality "
        "--S_algorithm skani "
        "--P_ani 0.9 "
        "--S_ani 0.95 "
        "--cov_thresh 0.1 "
        "--processors {threads} "
        "&> {log} "

#####################################
### Reference-based (GlobDB) runs ###
#####################################
# Same tools as above, but clustering sample genomes together with the full
# GlobDB r232 reference set in one run.
rule galah_globdb:
    input:
        "genome_lists_with_reference/genomes_{genome_count}_rep{repeat}.txt",
    output:
        directory("galah_globdb/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/galah_globdb/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/galah_globdb/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e galah "
        "galah cluster "
        "--genome-fasta-list {input} "
        "--ani 95 "
        "--min-aligned-fraction 15 "
        "--skip-sanitise-headers "
        "--output-cluster-definition {output}/clusters.tsv "
        "--threads {threads} "
        "&> {log} "

rule galah_low_memory_globdb:
    input:
        "genome_lists_with_reference/genomes_{genome_count}_rep{repeat}.txt",
    output:
        directory("galah_low_memory_globdb/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/galah_low_memory_globdb/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/galah_low_memory_globdb/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e galah "
        "galah cluster "
        "--genome-fasta-list {input} "
        "--ani 95 "
        "--min-aligned-fraction 15 "
        "--skip-sanitise-headers "
        "--output-cluster-definition {output}/clusters.tsv "
        "--threads {threads} "
        "--low-memory "
        "&> {log} "

rule skder_greedy_globdb:
    input:
        "genome_links/genomes_{genome_count}_rep{repeat}",
    output:
        directory("skder_greedy_globdb/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/skder_greedy_globdb/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/skder_greedy_globdb/cluster_{genome_count}_rep{repeat}.txt"
    params:
        globdb_dirs = " ".join(GLOBDB_CHUNK_DIRS),
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e skder "
        "skder "
        "--genomes {input} {params.globdb_dirs} "
        "--percent-identity-cutoff 95 "
        "--aligned-fraction-cutoff 15 "
        "--dereplication-mode greedy "
        "--symlink "
        "--output-directory {output}/skder "
        "--threads {threads} "
        "&> {log} "

rule drep_skani_globdb:
    input:
        "genome_lists_with_reference/genomes_{genome_count}_rep{repeat}.txt",
    output:
        directory("drep_skani_globdb/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/drep_skani_globdb/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/drep_skani_globdb/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e drep "
        "dRep dereplicate "
        "{output} "
        "--genomes {input} "
        "--ignoreGenomeQuality "
        "--S_algorithm skani "
        "--P_ani 0.9 "
        "--S_ani 0.95 "
        "--cov_thresh 0.1 "
        "--processors {threads} "
        "&> {log} "

# galah's dedicated reference-clustering mode: sample genomes are clustered
# only against the (pre-clustered) GlobDB reference set, rather than combining
# everything into one all-vs-all run.
rule galah_reference:
    input:
        genomes = "genome_lists/genomes_{genome_count}_rep{repeat}.txt",
        reference = "genome_lists/globdb_reference.txt",
    output:
        directory("galah_reference/cluster_{genome_count}_rep{repeat}"),
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "logs/galah_reference/logs/cluster_{genome_count}_rep{repeat}.log"
    benchmark:
        "benchmarks/galah_reference/cluster_{genome_count}_rep{repeat}.txt"
    shell:
        "mkdir -p {output} && "
        "TMPDIR={TMPDIR} "
        "pixi run -e galah "
        "galah cluster "
        "--genome-fasta-list {input.genomes} "
        "--reference-genomes-list {input.reference} "
        "--ani 95 "
        "--min-aligned-fraction 15 "
        "--skip-sanitise-headers "
        "--output-cluster-definition {output}/clusters.tsv "
        "--threads {threads} "
        "&> {log} "

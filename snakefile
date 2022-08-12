# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

#configfile: "config/config.yaml"

import os

onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')
    print(f"Env TMPDIR={os.environ.get('TMPDIR', '<n/a>')}")

# Define samples from data directory using wildcards
SAMPLES, = glob_wildcards('fastq/{samples}.fastq.gz')

# WC Sanity Check :D
print("Found:")
for WLDCRD in SAMPLES:
    print(WLDCRD)
print("")

rule target:
    input:
        expand('2_humann3RumFunc/{samples}/{samples}_kneaddata_genefamilies.tsv', samples=SAMPLES),
        expand('2_humann3Ovine/{samples}/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam_genefamilies.tsv', samples=SAMPLES)

rule kneaddata:
    input:
        reads = 'fastq/{samples}.fastq.gz',
    output:
        outDir = directory('1_kneaddata/{samples}'),
        clnReads = temp('1_kneaddata/{samples}/{samples}_kneaddata.fastq'),
        ovineReads = temp('1_kneaddata/{samples}/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam.fastq'),
        silvaReads = temp('1_kneaddata/{samples}/{samples}_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq'),
        trpReads = temp('1_kneaddata/{samples}/{samples}_kneaddata.repeats.removed.fastq'),
        trimReads = temp('1_kneaddata/{samples}/{samples}_kneaddata.trimmed.fastq'),
        readStats = '1_kneaddata/{samples}.read.stats.txt}'
    log:
        'logs/{samples}.kneaddata.log'
    conda:
        'biobakery'
    threads:12
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--input {input.reads} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '--trimmomatic /home/perrybe/conda-envs/biobakery/share/trimmomatic '
        '--sequencer-source TruSeq3 '
        '-db ref/ARS_UI_Ramb_v2 '
        '-db ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA '
        '-o {output.outDir} && '
        'seqkit stats -j 12 -a 1_kneaddata/{output.outDir}/*.fastq > {output.readStats} '

rule human3RumFunc:
    input:
        clnReads = '1_kneaddata/{samples}/{samples}_kneaddata.fastq',
    output:
        rumFuncDir = directory('2_humann3RumFunc/{samples}'),
        genes = '2_humann3RumFunc/{samples}/{samples}_kneaddata_genefamilies.tsv',
        pathways = '2_humann3RumFunc/{samples}/{samples}_kneaddata_pathabundance.tsv',
        pathwaysCoverage = '2_humann3RumFunc/{samples}/{samples}_kneaddata_pathcoverage.tsv'
    log:
        'logs/{samples}.human3.RumFunc.log'
    conda:
        'biobakery'
    threads:12
    message:
        'humann3 profiling RumFunc: {wildcards.samples}\n'
    shell:
        'humann '
        '--threads {threads} '
        '--input {input.clnReads} '
        '--output {output.rumFuncDir} '
        '--bypass-nucleotide-search '
        '--memory-use maximum '
        '--input-format fastq '
        '--search-mode uniref90 '
        '--verbose '
        '--log-level INFO '
        '--o-log {log} '
        '--remove-temp-output '

rule human3Ovine:
    input:
        ovineReads = '1_kneaddata/{samples}/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam.fastq',
    output:
        OvineDir = directory('2_humann3Ovine/{samples}'),
        genes = '2_humann3Ovine/{samples}/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam_genefamilies.tsv',
        pathways = '2_humann3Ovine/{samples}/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam_pathabundance.tsv',
        pathwaysCoverage = '2_humann3Ovine/{samples}/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam_pathcoverage.tsv'
    log:
        'logs/{samples}.human3.Ovine.log'
    conda:
        'biobakery'
    threads:12
    message:
        'humann3 profiling RumFunc: {wildcards.samples}\n'
    shell:
        'humann '
        '--threads {threads} '
        '--input {input.ovineReads} '
        '--output {output.OvineDir} '
        '--bypass-nucleotide-search '
        '--memory-use maximum '
        '--input-format fastq '
        '--search-mode uniref90 '
        '--verbose '
        '--log-level INFO '
        '--o-log {log} '
        '--remove-temp-output '


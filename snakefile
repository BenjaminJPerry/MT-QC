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
SAMPLES, = glob_wildcards('fastq/{samples}_L001_R1_001.fastq.gz')

rule target:
    input:
        expand('1_kneaddata/{samples}/seq_kneaddata_paired_1.fastq', samples=SAMPLES),
        expand('1_kneaddata/{samples}/seq_kneaddata_paired_2.fastq', samples=SAMPLES)


rule kneaddata:
    input:
        read1='fastq/{samples}_L001_R1_001.fastq.gz',
        read2='fastq/{samples}_L001_R2_001.fastq.gz',
        ovineDB='ref/ARS_UI_Ramb_v2',
        silvaDB='ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA'
    output:
        sampleDir=directory('1_kneaddata/{samples}'),
        clnRead1='1_kneaddata/{samples}/seq_kneaddata_paired_1.fastq',
        clnRead2='1_kneaddata/{samples}/seq_kneaddata_paired_2.fastq'
    log:
        'logs/{samples}.kneaddata.log'
    conda:
        'biobakery'
    threads:12
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--input1 {input.read1} '
        '--input2 {input.read2} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '-db {input.ovineDB} '
        '-db {input.silvaDB} '
        '-o {output.sampleDir} '
#
# rule porechop:
#     input:
#         'fastq/{samples}.fastq.gz'
#     output:
#         '01_chop/{samples}.chop.fastq.gz'
#     log:
#         'logs/{samples}/porechop.log'
#     conda:
#         'porechop'
#     params:
#         checks=config['porechop']['check_reads'],
#         adpthresh=config['porechop']['adapter_threshold'],
#         midthresh=config['porechop']['middle_threshold'],
#         splitlen=config['porechop']['min_split_size']
#     resources:
#         tempdir=config['TMPDIR']
#     threads:4
#     message:
#         'Chopping: {wildcards.samples}\n'
#         'TMPDIR: {resources.tempdir}'
#     shell:
#         'porechop '
#         '--threads {threads} '
#         '--verbosity 1 '
#         '--format fastq.gz '
#         '--check_reads {params.checks} '
#         '--adapter_threshold {params.adpthresh} '
#         '--middle_threshold {params.midthresh} '
#         '--min_split_read_size {params.splitlen} '
#         '-i {input} '
#         '-o {output} '
#         '2>&1 | tee {log}'
#
# rule cutadapt:
#     input:
#         rules.porechop.output
#     output:
#         '02_cutadapt/{samples}.chop.primer.fastq.gz'
#     threads:4
#     log:
#         'logs/{samples}/cutadapt.primers.log'
#     conda:
#         'cutadapt'
#     params:
#         fwdPrimer=config['cutadapt']['fwd'],
#         revPrimer=config['cutadapt']['rev'],
#     resources:
#         tempdir=config['TMPDIR']
#     message:
#         'removing primers: {wildcards.samples}\n'
#         'TMPDIR: {resources.tempdir}'
#     shell:
#         'cutadapt '
#         '--discard-untrimmed '
#         '--action=retain '
#         '-j {threads} '
#         '--error-rate 0.2 '
#         '-g {params.fwdPrimer} '
#         '-a {params.revPrimer} '
#         '-o {output} '
#         '{input} '
#         '2>&1 | tee {log}'

configfile: "config.yaml"

import pandas as pd
import os
from glob import glob

# make sure the tmp directory exists
os.makedirs(config['tmp_dir'], exist_ok=True)

print("Loading run data...")
runs = pd.read_csv(config['runs_file'], index_col='run_id')
print("Loaded %d runs." % len(runs.index))

rule all:
    input:
        expand("qc/fastqc/{run_id}_fastqc.html",
               run_id=runs.index.values),
        "qc/alleyoop/k562_summary.tsv",
        expand("qc/alleyoop/rates/{run_id}.fastq_slamdunk_mapped_filtered_overallrates.csv",
               run_id=runs.index.values),
        expand("qc/alleyoop/utrrates/{run_id}.fastq_slamdunk_mapped_filtered_mutationrates_utr.csv",
               run_id=runs.index.values),
        expand("qc/alleyoop/tcperreadpos/{run_id}.fastq_slamdunk_mapped_filtered_tcperreadpos.csv",
               run_id=runs.index.values),
        expand("qc/alleyoop/tcperutrpos/{run_id}.fastq_slamdunk_mapped_filtered_tcperutr.csv",
               run_id=runs.index.values)

def get_srr(wildcards):
    return runs.srr[wildcards.run_id]

rule download_fastq:
    output:
        fq="data/fastq/{run_id}.fastq.gz"
    params:
        srr=get_srr,
        tmp_dir=config['tmp_dir'],
        fq=lambda wcs: config['tmp_dir'] + "/" + wcs.run_id
    threads: 8
    resources:
        mem_mb=1000
    conda:
        "envs/samtools.yaml"
    shell:
        """
        fasterq-dump -e {threads} -o {params.fq}.fastq -t {params.tmp_dir} {params.srr}
        bgzip -@ {threads} -c {params.fq}.fastq > {output.fq}
        rm {params.fq}.fastq
        """

rule fastqc:
    input:
        "data/fastq/{run_id}.fastq.gz"
    output:
        "qc/fastqc/{run_id}_fastqc.html"
    resources:
        mem_mb=1000
    conda:
        "envs/samtools.yaml"
    shell:
        """
        outdir=$(dirname {output})
        fastqc -o $outdir --nogroup {input}
        """

rule utrome_bed:
    input:
        gtf=config['utrome_gtf'],
        awk="scripts/utrome_gtf_to_bed.awk"
    output:
        "data/bed/utrome.bed"
    shell:
        """
        awk -f {input.awk} {input.gtf} > {output}
        """

rule slamdunk_all:
    input:
        fqs=expand("data/fastq/{run_id}.fastq.gz",
                   run_id=runs.index.values),
        csv=config['samples_file'],
        bed="data/bed/utrome.bed",
        fa=config['genome_fa']
    output:
        directory("data/slamdunk/count")
    params:
        rl=config['max_read_length']
    threads: 24
    resources:
        mem_mb=2000
    conda: "envs/slamdunk.yaml"
    shell:
        """
        slamdunk all \\
          -r {input.fa} \\
          -b {input.bed} \\
          -o data/slamdunk \\
          -5 12 -rl {params.rl} \\
          -n 100 -m \\
          -t {threads} \\
          --skip-sam \\
          {input.csv}
        """

rule alleyoop_summary:
    input:
        bams=expand("data/slamdunk/filter/{run_id}.fastq_slamdunk_mapped_filtered.bam",
                    run_id=runs.index.values),
        count_dir=rules.slamdunk_all.output
    output:
        tsv="qc/alleyoop/k562_summary.tsv",
        pdf="qc/alleyoop/k562_summary_PCA.pdf",
        pca="qc/alleyoop/k562_summary_PCA.txt"
    conda: "envs/slamdunk.yaml"
    shell:
        """
        alleyoop summary \\
          -o {output.tsv} \\
          -t {input.count_dir} \\
          {input.bams}
        """

rule alleyoop_rates:
    input:
        bam="data/slamdunk/filter/{run_id}.{filtered}.bam",
        fa=config['genome_fa']
    output:
        csv="qc/alleyoop/rates/{run_id}.{filtered}_overallrates.csv",
        pdf="qc/alleyoop/rates/{run_id}.{filtered}_overallrates.pdf"
    conda: "envs/slamdunk.yaml"
    wildcard_constraints:
        filtered='fastq_slamdunk_mapped_filtered'
    threads: 4
    resources:
        mem_mb=2000
    shell:
        """
        alleyoop rates \\
          -o qc/alleyoop/rates \\
          -r {input.fa} \\
          -t {threads} \\
          {input.bam}
        """

rule alleyoop_utrrates:
    input:
        bam="data/slamdunk/filter/{run_id}.{filtered}.bam",
        bed="data/bed/utrome.bed",
        fa=config['genome_fa']
    output:
        csv="qc/alleyoop/utrrates/{run_id}.{filtered}_mutationrates_utr.csv",
        pdf="qc/alleyoop/utrrates/{run_id}.{filtered}_mutationrates_utr.pdf"
    params:
        rl=config['max_read_length']
    conda: "envs/slamdunk.yaml"
    wildcard_constraints:
        filtered='fastq_slamdunk_mapped_filtered'
    threads: 4
    resources:
        mem_mb=2000
    shell:
        """
        alleyoop utrrates \\
          -o qc/alleyoop/utrrates \\
          -r {input.fa} \\
          -b {input.bed} \\
          -l {params.rl} \\
          -t {threads} \\
          {input.bam}
        """

rule alleyoop_tcperreadpos:
    input:
        bam="data/slamdunk/filter/{run_id}.{filtered}.bam",
        fa=config['genome_fa'],
        snp=directory("data/slamdunk/snp")
    output:
        csv="qc/alleyoop/tcperreadpos/{run_id}.{filtered}_tcperreadpos.csv",
        pdf="qc/alleyoop/tcperreadpos/{run_id}.{filtered}_tcperreadpos.pdf"
    params:
        rl=config['max_read_length']
    conda: "envs/slamdunk.yaml"
    wildcard_constraints:
        filtered='fastq_slamdunk_mapped_filtered'
    threads: 4
    resources:
        mem_mb=2000
    shell:
        """
        alleyoop tcperreadpos \\
          -o qc/alleyoop/tcperreadpos \\
          -r {input.fa} \\
          -l {params.rl} \\
          -s {input.snp} \\
          -t {threads} \\
          {input.bam}
        """

rule alleyoop_tcperutrpos:
    input:
        bam="data/slamdunk/filter/{run_id}.{filtered}.bam",
        fa=config['genome_fa'],
        bed="data/bed/utrome.bed",
        snp=directory("data/slamdunk/snp")
    output:
        csv="qc/alleyoop/tcperutrpos/{run_id}.{filtered}_tcperutr.csv",
        pdf="qc/alleyoop/tcperutrpos/{run_id}.{filtered}_tcperutr.pdf"
    params:
        rl=config['max_read_length']
    conda: "envs/slamdunk.yaml"
    wildcard_constraints:
        filtered='fastq_slamdunk_mapped_filtered'
    threads: 2
    resources:
        mem_mb=4000
    shell:
        """
        alleyoop tcperutrpos \\
          -o qc/alleyoop/tcperutrpos \\
          -r {input.fa} \\
          -b {input.bed} \\
          -l {params.rl} \\
          -s {input.snp} \\
          -t {threads} \\
          {input.bam}
        """

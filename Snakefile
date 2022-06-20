"""
Command to run:
 snakemake -pr -j 16 --use-conda --conda-frontend mamba --keep-going \
 --latency-wait 340 --rerun-incomplete --conda-prefix ~/envs/ \
 --cluster  "sbatch --time={resources.time} --mem={resources.mem} --ntasks={resources.threads} \
 --output log/%x.%j.out --partition=CLUSTER"

Where
mem: 80G (kaiju 250)
time: 12:00:00 h
ntasks: 8
"""

import os

localrules: mviest_summary, mviest_plot, check_input, count_virome_contigs, 
            filter_input_by_length, contig_summary, contig_selector, vs2_summary

configfile: "config.yaml"
### all ########################################################################
rule all:
    input:
        expand("results/{sample}/log/00_INPUT_CHECK_DONE.log", sample = config["samples"]),
        expand("results/{sample}/{sample}_mviest_summary.tsv", sample = config["samples"]),
        expand("results/{sample}/{sample}_mviest_plot.png", sample = config["samples"])


### summary steps ##############################################################
rule mviest_summary:
    input: 
        "results/{sample}/{sample}_mviest_plot.png"
    output:
        "results/{sample}/{sample}_mviest_summary.tsv"
    shell:
        "touch {output}"

rule mviest_plot:
    conda:
        "envs/R.yaml"
    input:
        kaiju_mvome_table="results/{sample}/kaiju/mvome_summary.tsv",
        kaiju_virome_table="results/{sample}/kaiju/virome_summary.tsv",
        viromeQC_mvome="results/{sample}/viromeQC/mvome_QC.tsv",
        viromeQC_virome="results/{sample}/viromeQC/virome_QC.tsv",
        contig_summary="results/{sample}/contig_summary_{sample}.tsv",
        mvome_reads_vs_metagenome_scafstats="results/{sample}/mappings/mvome_positive_vs_metagenome/scafstats.txt"
    output:
        "results/{sample}/{sample}_mviest_plot.png"
    shell:
        """
        Rscript scripts/mviest_plot.R \
        {input.kaiju_mvome_table} {input.kaiju_virome_table} \
        {input.viromeQC_mvome} {input.viromeQC_virome} \
        {input.contig_summary} \
        {input.mvome_reads_vs_metagenome_scafstats} \
        {output}
        """


### sanity checks ##############################################################
rule check_input:
    input:
        vr1="results/{sample}/reads/virome.reads1.fastq.gz",
        vr2="results/{sample}/reads/virome.reads2.fastq.gz",
        vc="results/{sample}/contigs/virome.contigs.fasta",
        mc="results/{sample}/contigs/metagenome.contigs.fasta"
    output:
        "results/{sample}/log/00_INPUT_CHECK_DONE.log"
    shell:
        "touch {output}"

rule count_virome_contigs:
    input:
        "results/{sample}/contigs/virome.contigs.filtered.fasta"
    output:
        enough="results/{sample}/enough_filtered_virome_contigs.txt"
    params:
        minimum_number_of_viral_contigs=config["minimum_number_of_virome_input_contigs"]
    shell:
        """
        n=$(grep '>' {input}  | wc -l)
        if [[ $n -ge {params.minimum_number_of_viral_contigs} ]]
        then
            touch {output.enough}
        else
            touch results/{wildcards.sample}/NOT_enough_filtered_virome_contigs.txt
        fi
        """


### handle contigs #############################################################
rule filter_input_virome_by_length:
    conda:
        "envs/python.yaml"
    input:
        vc="results/{sample}/contigs/virome.contigs.fasta",
        input_check="results/{sample}/log/00_INPUT_CHECK_DONE.log"
    output:
        vc_filtered="results/{sample}/contigs/virome.contigs.filtered.fasta",
    params:
        virome_min_len=config["virome_min_len"],
    shell:
        """
        python3 scripts/filter_by_len.py {input.vc} \
        {params.virome_min_len} > {output.vc_filtered}
        """

rule filter_input_metagenome_by_length:
    conda:
        "envs/python.yaml"
    input:
        mc="results/{sample}/contigs/metagenome.contigs.fasta",
        input_check="results/{sample}/log/00_INPUT_CHECK_DONE.log"
    output:
        mc_filtered="results/{sample}/contigs/metagenome.contigs.filtered.fasta"
    params:
        metagenome_min_min_len=config["metagenome_min_min_len"]
    shell:
        """
        python3 scripts/filter_by_len.py {input.mc} \
        {params.metagenome_min_min_len} > {output.mc_filtered}
        """

rule contig_summary:
    conda:
        "envs/R.yaml"
    input:
        contigs="results/{sample}/contigs/virome.contigs.filtered.fasta",
        vs2_out="results/{sample}/vs2_final_report.tsv",
        dvf_out="results/{sample}/dvf/virome.contigs.filtered.fasta_gt1bp_dvfpred.txt",
        platon_out="results/{sample}/platon_out/platon_plasmids.log",
        enough="results/{sample}/enough_filtered_virome_contigs.txt"
    output:
        "results/{sample}/contig_summary_{sample}.tsv"
    params:
        cutoff=config["long_contig_cutoff"]
    log:
        "results/{sample}/log/text_summary.log"
    shell:
        "Rscript scripts/summarize.R {wildcards.sample} {params.cutoff} 2> {log}"

rule contig_selector:
    conda:
        "envs/python.yaml"
    input:
        "results/{sample}/contig_summary_{sample}.tsv",
        "results/{sample}/contigs/virome.contigs.filtered.fasta",
        "results/{sample}/enough_filtered_virome_contigs.txt"
    output:
        mv_contigs="results/{sample}/contigs/true.mvome.fasta",
        virome_contigs="results/{sample}/contigs/true.virome.fasta"
    log:
        "results/{sample}/log/contig_selector.log"
    shell:
        "python scripts/contig_selector.py {wildcards.sample} 2> {log}"


### kaiju ######################################################################
rule kaiju:
    conda:
        "envs/kaiju.yaml"
    input:
        mvome_r1="results/{sample}/reads/true.mvome.reads1.fq.gz",
        mvome_r2="results/{sample}/reads/true.mvome.reads2.fq.gz",
        virome_r1="results/{sample}/reads/true.virome.reads1.fq.gz",
        virome_r2="results/{sample}/reads/true.virome.reads2.fq.gz"
    output:
        mvome_raw="results/{sample}/kaiju/mvome.out",
        virome_raw="results/{sample}/kaiju/virome.out",
        mvome_table="results/{sample}/kaiju/mvome_summary.tsv",
        virome_table="results/{sample}/kaiju/virome_summary.tsv"
    params:
        kaiju_fmi=config["kaiju_fmi"],
        kaiju_nodes=config["kaiju_nodes"],
        kaiju_names=config["kaiju_names"]
    resources:
        mem=config["sbatch_high_mem"],
        time="12:00:00",
        ntasks=8
    shell:
        """
        mkdir -p results/{wildcards.sample}/kaiju
        
        kaiju-multi -z {resources.ntasks} \
        -t {params.kaiju_nodes} -f {params.kaiju_fmi} \
        -i {input.mvome_r1},{input.virome_r1} \
        -j {input.mvome_r2},{input.virome_r2} \
        -o {output.mvome_raw},{output.virome_raw}

        kaiju2table -t {params.kaiju_nodes} -n {params.kaiju_names} \
        -r genus -o {output.mvome_table} {output.mvome_raw}

        kaiju2table -t {params.kaiju_nodes} -n {params.kaiju_names} \
        -r genus -o {output.virome_table} {output.virome_raw}

        """


### viromeQC ###################################################################
rule viromeQC:
    """
    This runs viromeQC on the input reads, mvome reads and creates a report file.
    """
    conda:
        "envs/viromeQC.yaml"
    params:
        viromeQC_path=config["viromeQC_path"]
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    input:
        mvome_r1="results/{sample}/reads/true.mvome.reads1.fq.gz",
        mvome_r2="results/{sample}/reads/true.mvome.reads2.fq.gz",
        virome_r1="results/{sample}/reads/true.virome.reads1.fq.gz",
        virome_r2="results/{sample}/reads/true.virome.reads2.fq.gz"
    output:
        mvome_out="results/{sample}/viromeQC/mvome_QC.tsv",
        virome_out="results/{sample}/viromeQC/virome_QC.tsv"        
    shell:
        """
        mkdir -p results/{wildcards.sample}/viromeQC

        python {params.viromeQC_path}/viromeQC.py -w environmental \
        --diamond_threads {resources.ntasks} --bowtie2_threads {resources.ntasks} \
        -i {input.mvome_r1} {input.mvome_r2} \
        -o {output.mvome_out}

        python {params.viromeQC_path}/viromeQC.py -w environmental \
        --diamond_threads {resources.ntasks} --bowtie2_threads {resources.ntasks} \
        -i {input.virome_r1} {input.virome_r2} \
        -o {output.virome_out}
        """


### viral prediction ###########################################################
rule DeepVirFinder:
    conda:
        "envs/dvf.yaml"
    input: 
        contigs="results/{sample}/contigs/virome.contigs.filtered.fasta",
        enough="results/{sample}/enough_filtered_virome_contigs.txt"
    output:
        "results/{sample}/dvf/virome.contigs.filtered.fasta_gt1bp_dvfpred.txt"
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    params:
        dvf_executable=config["dvf_path"]
    shell:
        """
        python {params.dvf_executable} \
        -i results/{wildcards.sample}/contigs/virome.contigs.filtered.fasta \
        -o results/{wildcards.sample}/dvf -c {resources.ntasks}
        """

rule vs2:
    conda:
        "envs/vs2.yaml"
    input:
        virome_contigs="results/{sample}/contigs/virome.contigs.filtered.fasta",
        enough="results/{sample}/enough_filtered_virome_contigs.txt"
    output:
        "results/{sample}/vs2/final-viral-combined.fa",
        "results/{sample}/vs2/final-viral-score.tsv"
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    params:
        vsdbdir=config["vs2_db"]
    shell:
        """
        virsorter run --keep-original-seq -i {input.virome_contigs} \
        -w results/{wildcards.sample}/vs2 \
        --min-score 0.5 -j {resources.ntasks} all --db-dir {params.vsdbdir}
        """

rule checkv:
    conda:
        "envs/vs2.yaml"
    input:
        "results/{sample}/vs2/final-viral-combined.fa"
    output:
        "results/{sample}/checkv/contamination.tsv",
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    params:
        checkvdbdir=config["checkv_db"]
    shell:
        """
        checkv end_to_end {input} results/{wildcards.sample}/checkv \
        -t {resources.ntasks} -d {params.checkvdbdir}
        """

rule vs2_summary:
    conda:
        "envs/R.yaml"
    input:
        "results/{sample}/vs2/final-viral-score.tsv",
        "results/{sample}/checkv/contamination.tsv"
    output:
        "results/{sample}/vs2_final_report.tsv"
    shell:
        "Rscript scripts/analyze_vs2_output.R {input} {output}"


### plasmid prediction #########################################################
rule platon:
    conda:
        "envs/platon.yaml"
    input:
        virome_contigs="results/{sample}/contigs/virome.contigs.filtered.fasta",
        enough="results/{sample}/enough_filtered_virome_contigs.txt"
    output:
        "results/{sample}/platon_out/platon_plasmids.log"
    params:
        platon_db=config["platon_db"]
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    shell:
        """
        platon --db {params.platon_db} --prefix platon_plasmids \
        --output results/{wildcards.sample}/platon_out/ \
        --threads {resources.ntasks} \
        {input.virome_contigs}
        """


### mapping ####################################################################
rule map_virome_reads_to_mv_contigs:
    conda:
        "envs/bbmap.yaml" # "echo lala"
    input:
        r1="results/{sample}/reads/virome.reads1.fastq.gz",
        r2="results/{sample}/reads/virome.reads2.fastq.gz",
        mv_contigs="results/{sample}/contigs/true.mvome.fasta"
    output:
        r1="results/{sample}/reads/true.mvome.reads1.fq.gz",
        r2="results/{sample}/reads/true.mvome.reads2.fq.gz",
        statsfile="results/{sample}/mappings/input_reads_vs_mv_contigs/statsfile.txt",
    log:
        "results/{sample}/log/map_virome_reads_to_mv_contigs.log"
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    shell:
        """
        mkdir -p results/{wildcards.sample}/mappings/input_reads_vs_mv_contigs/

        bbmap.sh in={input.r1} in2={input.r2} \
        ref={input.mv_contigs} \
        outm={output.r1} outm2={output.r2} \
        statsfile={output.statsfile} \
        t={resources.ntasks} \
        path=results/{wildcards.sample}/
        """

rule map_virome_reads_to_true_virome_contigs:
    conda:
        "envs/bbmap.yaml" # "echo lala"
    input:
        r1="results/{sample}/reads/virome.reads1.fastq.gz",
        r2="results/{sample}/reads/virome.reads2.fastq.gz",
        virome_contigs="results/{sample}/contigs/true.virome.fasta"
    output:
        r1="results/{sample}/reads/true.virome.reads1.fq.gz",
        r2="results/{sample}/reads/true.virome.reads2.fq.gz",
        scafstats="results/{sample}/mappings/input_reads_vs_true_virome_contigs/scafstats.txt",
        rpkm="results/{sample}/mappings/input_reads_vs_true_virome_contigs/rpkm.txt"
    log:
        "results/{sample}/log/map_virome_reads_vs_true_virome_contigs.log"
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    shell:
        """
        mkdir -p results/{wildcards.sample}/mappings/input_reads_vs_true_virome_contigs/

        bbmap.sh in={input.r1} in2={input.r2} \
        ref={input.virome_contigs} \
        outm={output.r1} outm2={output.r2} \
        scafstats={output.scafstats} \
        rpkm={output.rpkm} \
        t={resources.ntasks} \
        path=results/{wildcards.sample}/
        """

rule map_mv_positive_reads_to_metagenome:
    conda:
        "envs/bbmap.yaml" # "echo lala"
    input:
        r1="results/{sample}/reads/true.mvome.reads1.fq.gz",
        r2="results/{sample}/reads/true.mvome.reads2.fq.gz",
        metagenome_contigs="results/{sample}/contigs/metagenome.contigs.filtered.fasta"
    output:
        rpkm="results/{sample}/mappings/mvome_positive_vs_metagenome/rpkm.txt",
        scafstats="results/{sample}/mappings/mvome_positive_vs_metagenome/scafstats.txt"
    resources:
        mem=config["sbatch_normal_mem"],
        time="12:00:00",
        ntasks=8
    shell:
        """
        mkdir -p results/{wildcards.sample}/mappings/mvome_positive_vs_metagenome/

        bbmap.sh in={input.r1} in2={input.r2} \
        ref={input.metagenome_contigs} \
        scafstats={output.scafstats} \
        rpkm={output.rpkm} \
        t={resources.ntasks} \
        path=results/{wildcards.sample}/
        """    
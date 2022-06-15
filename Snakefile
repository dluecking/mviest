import os

localrules: filter_input_by_length, count_virome_contigs,
            check_input, vs2_summary

configfile: "config.yaml"
### all ########################################################################
rule all:
    input:
        expand("results/{sample}/{sample}_mvip_summary.tsv", sample = config["samples"]),
        expand("results/{sample}/{sample}_mvip_plot.png", sample = config["samples"])



### summary steps ##############################################################
rule mvip_summary:
    input: 
    
    output:
        "results/{sample}/{sample}_mvip_summary.tsv"
    shell:

rule mvip_plot:
    input: 
    
    output:
        "results/{sample}/{sample}_mvip_plot.png""
    shell:



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
rule filter_input_by_length:
    conda:
        "envs/python.yaml"
    input:
        vc="results/{sample}/contigs/virome.contigs.fasta",
        mc="results/{sample}/contigs/metagenome.contigs.fasta",
        input_check="results/{sample}/log/00_INPUT_CHECK_DONE.log"
    output:
        vc_filtered="results/{sample}/contigs/virome.contigs.filtered.fasta",
        mc_filtered="results/{sample}/contigs/metagenome.contigs.filtered.fasta"
    params:
        virome_min_len=config["virome_min_len"],
        metagenome_min_min_len=config["metagenome_min_min_len"]
    shell:
        """
        python3 scripts/filter_by_len.py {input.vc} \
        {params.virome_min_len} > {output.vc_filtered}

        python3 scripts/filter_by_len.py {input.mc} \
        {params.virome_min_len} > {output.mc_filtered}
        """

rule contig_summary:
    conda:
        "envs/R.yaml"
    input:
        contigs="results/{sample}/contigs/virome.contigs.filtered.fasta",
        vs2_out="results/{sample}/vs2_final_report.tsv",
        dvf_out="results/{sample}/dvf/virome.contigs.filtered.fasta_gt2000bp_dvfpred.txt",
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


### viromeQC ###################################################################
rule viromeQC:
    """
    This runs viromeQC on the input reads, mvome reads and creates a report file.
    """
    conda:
        "envs/viromeQC.yaml"
    params:
        viromeQC_path=config["viromeQC_path"]
    threads: 8
    input:
        mvome_r1="results/{sample}/reads/mvome.reads1.fq.gz",
        mvome_r2="results/{sample}/reads/mvome.reads2.fq.gz",
        virome_r1="results/{sample}/reads/mvome.reads2.fq.gz",
        virome_r2="results/{sample}/reads/mvome.reads2.fq.gz"
    output:
        mvome_out="results/{sample}/viromeQC_mvome_{sample}.tsv",
        virome_out="results/{sample}/viromeQC_virome_{sample}.tsv"        
    shell:
        """
        python {params.viromeQC_path}/viromeQC.py -w environmental \
        --diamond_threads {threads} --bowtie2_threads {threads} \
        -i {input.mvome_r1} {input.mvome_r2} \
        -o {output.mvome_out}

        python {params.viromeQC_path}/viromeQC.py -w environmental \
        --diamond_threads {threads} --bowtie2_threads {threads} \
        -i {input.virome_r1} {input.virome_r2} \
        -o {output.virome_out}
        """




### viral prediction ###########################################################
rule DeepVirFinder:
    conda:
        "envs/dvf.yaml"
    input: 
        "results/{sample}/contigs/virome.contigs.filtered.fasta"
    output:
        "results/{sample}/dvf/virome.contigs.filtered.fasta_gt2000bp_dvfpred.txt"
    threads: 8
    params:
        dvf_executable=config["dvf_path"]
    shell:
        """
        python {params.dvf_executable} \
        -i results/{wildcards.sample}/contigs/virome.contigs.filtered.fasta \
        -o results/{wildcards.sample}/dvf -c {threads}
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
    threads: 8
    params:
        vsdbdir=config["vs2_db"]
    shell:
        """
        virsorter run --keep-original-seq -i {input.virome_contigs} \
        -w results/{wildcards.sample}/vs2-pass1 \
        --min-score 0.5 -j {threads} all --db-dir {params.vsdbdir}
        """

rule checkv:
    conda:
        "envs/vs2.yaml"
    input:
        "results/{sample}/vs2/final-viral-combined.fa"
    output:
        "results/{sample}/checkv/contamination.tsv",
    threads: 8
    params:
        checkvdbdir=config["checkv_db"]
    shell:
        """
        checkv end_to_end {input} results/{wildcards.sample}/checkv \
        -t {threads} -d {params.checkvdbdir}
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
    threads: 8
    shell:
        """
        platon --db {params.platon_db} --prefix platon_plasmids \
        --output results/{wildcards.sample}/platon_out/ \
        --threads {threads}
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
        r1="results/{sample}/reads/mvome.reads1.fq.gz",
        r2="results/{sample}/reads/mvome.reads2.fq.gz",
        statsfile="results/{sample}/mappings/input_reads_vs_mv_contigs/statsfile.txt",
    log:
        "results/{sample}/log/map_virome_reads_to_mv_contigs.log"
    threads: 8
    shell:
        """
        mkdir -p results/{wildcards.sample}/mappings/input_reads_vs_mv_contigs/

        bbmap.sh in={input.r1} in2={input.r2} \
        ref={input.mv_contigs} \
        outm={output.r1} outm2={output.r2} \
        statsfile={output.statsfile} \
        """

rule map_virome_reads_to_true_virome_contigs:
    conda:
        "envs/bbmap.yaml" # "echo lala"
    input:
        r1="results/{sample}/reads/virome.reads1.fastq.gz",
        r2="results/{sample}/reads/virome.reads2.fastq.gz",
        virome_contigs="results/{sample}/contigs/true.virome.fasta"
    output:
        scafstats="results/{sample}/mappings/input_reads_vs_true_virome_contigs/scafstats.txt",
        rpkm="results/{sample}/mappings/input_reads_vs_true_virome_contigs/rpkm.txt"
    log:
        "results/{sample}/log/map_virome_reads_vs_true_virome_contigs.log"
    threads: 8
    shell:
        """
        mkdir -p results/{wildcards.sample}/mappings/input_reads_vs_true_virome_contigs/

        bbmap.sh in={input.r1} in2={input.r2} \
        ref={input.virome_contigs} \
        scafstats={output.statsfile} \
        """

rule map_mv_positive_reads_to_mv_producer:
    conda:
        "envs/bbmap.yaml" # "echo lala"
    params:
        reference_list = config["mv_producer_list"]
    input:
        r1="results/{sample}/reads/mvome.reads1.fq.gz",
        r2="results/{sample}/reads/mvome.reads2.fq.gz"
    output:
        "results/{sample}/all_mv_producer_mapped.txt"
    threads: 8
    shell:
        """
        mkdir -p results/{wildcards.sample}/mappings/mv_producer

        for reference in `cat {params.reference_list}`; do
            bash scripts/map_reads_to_reference.sh \
            {input.r1} data/mv_producer/genomes/${{reference}} \
            {wildcards.sample}_r1_vs_${{reference}} \
            results/{wildcards.sample}/mappings/mv_producer {threads};

            bash scripts/map_reads_to_reference.sh \
            {input.r2} data/mv_producer/genomes/${{reference}} \
            {wildcards.sample}_r2_vs_${{reference}} \
            results/{wildcards.sample}/mappings/mv_producer {threads};
        done;

        rm results/{wildcards.sample}/mappings/mv_producer/*.fastq.gz
        rm results/{wildcards.sample}/mappings/mv_producer/*.bam
        touch {output}
        """    

rule chromap_index:
    """
    Build chromap index
    """
    input:
        "genomes/genome.fasta.gz"
    output:
        "genomes/genome.index"
    log:
        chromap = "logs/chromap_idx.log",
    threads:
        max_threads
    resources:
        mem_mb = 50000
    conda:
        "../envs/chromap.yaml"
    shell:
        "chromap -i -r {input} -o {output} 2> "
        "{log.chromap}"

rule chromap:
    """
    Run Chromap
    """
    input:
        fastq_1 = "fastqs/{sample}/r1.fastq.gz",
        fastq_2 = "fastqs/{sample}/r2.fastq.gz",
        fastq_bc = "fastqs/{sample}/bc.fastq.gz",
        wl = "bc_whitelist.txt",
        ref = "genomes/genome.fasta.gz",
        index = "genomes/genome.index"
    output:
        "results/{sample}/alignments_unsorted.sam"
    log:
        chromap = "logs/{sample}/chromap.log",
    threads:
        max_threads
    resources:
        runtime_min = 480
    conda:
        "../envs/chromap.yaml"
    shell:
        "chromap --preset atac --SAM --drop-repetitive-reads 4 -q 0 --trim-adapters -x {input.index} -t {threads} "
        "-r {input.ref} -1 {input.fastq_1} -2 {input.fastq_2} -o {output} -b {input.fastq_bc} --barcode-whitelist {input.wl} 2> "
        "{log.chromap}"

rule collate_alignments:
    """
    Run fixmate on alignments
    """
    input:
        "results/{sample}/alignments_unsorted.sam"
    output:
        "results/{sample}/alignments_collated.bam"
    threads:
        max_threads
    resources:
        runtime_min = 480
    conda:
        "../envs/chromap.yaml"
    shell:
        "samtools collate -@ {threads} -o {output} {input}"

rule fixmate:
    """
    Run fixmate on alignments
    """
    input:
        "results/{sample}/alignments_collated.bam"
    output:
        "results/{sample}/alignments_fixmate.bam"
    resources:
        mem_mb = 1000
    conda:
        "../envs/chromap.yaml"
    shell:
        "samtools fixmate -r {input} {output}"

rule filter_mito:
    """
    Filter and count mitochondrial reads (and also fiter out secondary alignments)
    """
    input: 
        "results/{sample}/alignments_fixmate.bam"
    output: 
        bam = "results/{sample}/alignments_no_mito.bam",
        qc = "results/{sample}/frac_mito.tsv"
    params:
        mitochr = lambda w: config["mito_chr"]
    resources:
        mem_mb = 1000
    conda:
        "../envs/chromap.yaml"
    script:
        "../scripts/filter_mito.py"

rule sort_alignments:
    """
    Sort alignments
    """
    input:
        "results/{sample}/alignments_no_mito.bam"
    output:
        "results/{sample}/alignments_sorted.bam"
    threads:
        max_threads
    resources:
        runtime_min = 480
    conda:
        "../envs/chromap.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"


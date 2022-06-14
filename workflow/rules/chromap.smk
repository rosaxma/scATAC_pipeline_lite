
rule chromap_index:
    """
    Build chromap index
    """
    input:
        "genomes/genome.fa"
    output:
        "genomes/genome.index"
    log:
        chromap = "logs/chromap_idx.log",
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
        fastq_1 = "fastqs/{sample}/chromap/r1.fastq.gz",
        fastq_2 = "fastqs/{sample}/chromap/r2.fastq.gz",
        fastq_bc = "fastqs/{sample}/chromap/bc.fastq.gz",
        wl = "bc_whitelist.txt",
        ref = "genomes/genome.fasta.gz",
        index = "genomes/genome.index"
    output:
        "results/{sample}/chromap/alignments_unsorted.bam"
    params:
        barcode_dist = lambda w: config["max_barcode_dist"],
        multimapping = config["multimapping"]
    log:
        chromap = "logs/{sample}/chromap/chromap.log"
    threads:
        max_threads
    resources:
        runtime_min = 480
    conda:
        "../envs/chromap.yaml"
    shell:
        "chromap --preset atac --SAM --drop-repetitive-reads {params.multimapping} -q 0 --trim-adapters -t {threads} --bc-error-threshold {params.barcode_dist} "
        "-x {input.index} -r {input.ref} -1 {input.fastq_1} -2 {input.fastq_2} -o /dev/stdout -b {input.fastq_bc} --barcode-whitelist {input.wl} 2> "
        "{log.chromap} | "
        "samtools view -b -S -o {output} -"

rule collate_alignments:
    """
    Run fixmate on alignments
    """
    input:
        "results/{sample}/chromap/alignments_unsorted.bam"
    output:
        "results/{sample}/chromap/alignments_collated.bam"
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
        "results/{sample}/chromap/alignments_collated.bam"
    output:
        "results/{sample}/chromap/alignments_map_out.bam"
    resources:
        mem_mb = 1000
    conda:
        "../envs/chromap.yaml"
    shell:
        "samtools fixmate -r {input} {output}"






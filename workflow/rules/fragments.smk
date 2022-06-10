rule bam_to_fragments: 
    """
    Convert BAM to fragment file
    """
    input:
        "results/{sample}/alignments_sorted.bam"
    output:
        "results/{sample}/fragments.tsv"
    params:
        shift_plus = config["tn5_shift_plus"],
        shift_minus = config["tn5_shift_minus"]
    resources:
        mem_mb = 1000
    conda:
        "../envs/fragments.yaml"
    script:
        "../scripts/bam_to_fragments.py"

rule compress_fragments:
    """
    Compress fragment file
    """
    input:
        "results/{sample}/fragments.tsv"
    output: 
        "results/{sample}/fragments.tsv.gz"
    resources:
        mem_mb = 1000
    conda:
        "../envs/fragments.yaml"
    shell: 
        "bgzip -c {input} > {output}"

rule index_fragments:
    """
    Index fragment file
    """
    input: 
        "results/{sample}/fragments.tsv.gz"
    output: 
        "results/{sample}/fragments.tsv.gz.tbi"
    resources:
        mem_mb = 1000
    conda:
        "../envs/fragments.yaml"
    shell: 
        "tabix --zero-based --preset bed {input}"
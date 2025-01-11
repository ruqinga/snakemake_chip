rule bam2bw:
    input:
        bam = get_uniq_list
    output:
        bw="{bw_out}/{sample}_unique_{binsize}.bw"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        bw_out = directories["bw_out"],
        thread = config["threads"]
    shell:
        """
        # bam2bw
        bamCoverage --bam {input.bam} --outFileName {output.bw} --outFileFormat bigwig \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread}

        echo "visualization finished at $(date)"
        """

rule bam2bed:
    input:
        bam = get_uniq_list
    output:
        bed="{bw_out}/{sample}_unique_{binsize}.bed"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        bw_out = directories["bw_out"],
        thread = config["threads"]
    shell:
        """
        # bam2bw
        bamCoverage --bam {input.bam} --outFileName {output.bed} --outFileFormat bedgraph \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread}

        echo "visualization finished at $(date)"
        """
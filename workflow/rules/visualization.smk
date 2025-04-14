rule bam2bw:
    input:
        bam = get_uniq_list
    output:
        bw="Results/06_visualization/{sample}_unique_{binsize}.bw"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        thread = config["threads"]
    log:
        log = "Results/06_visualization/logs/{sample}_{binsize}.log"
    shell:
        """
        # bam2bw
        bamCoverage --bam {input.bam} --outFileName {output.bw} --outFileFormat bigwig \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread} >> {log.log} 2>&1
        """

rule bam2bed:
    input:
        bam = get_uniq_list
    output:
        bed="Results/06_visualization/bed/{sample}_unique_{binsize}.bed"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        thread = config["threads"]
    log:
        log = "Results/06_visualization/logs/{sample}_{binsize}.log"
    shell:
        """
        # bam2bed
        bamCoverage --bam {input.bam} --outFileName {output.bed} --outFileFormat bedgraph \
              --binSize {wildcards.binsize} --normalizeUsing RPKM --numberOfProcessors {params.thread} >> {log.log} 2>&1
        """
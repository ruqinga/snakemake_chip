rule bam2bw:
    input:
        uniq_bam = get_uniq_list
    output:
        bw = "{align_out}/{sample}.sam",
        sorted_bam = "{align_out}/{sample}.Hisat_aln.sorted.bam",
        sorted_bam_bai = "{align_out}/{sample}.Hisat_aln.sorted.bam.bai"
    conda:
        config["conda_env"]
    group: "processing"
    priority: 2
    params:
        option = config["hisat2"]["params"],
        index= config["hisat2"]["index"][config["Spe"]],
        align_out = directories["align_out"],
        log_dir= directories["log_out"]
    shell:
        """
        bamCoverage --bam "{input.uniq_bam}" --outFileName ${workdir}/bam2bw/${prefix}.bw --outFileFormat bigwig \
              --binSize 100 --normalizeUsing RPKM --numberOfProcessors 20 >${workdir}/bam2bw/${prefix}.log 2>&1
        """

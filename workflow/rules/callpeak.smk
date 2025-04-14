rule callpeak_pe:
    input:
        bam = get_uniq_list
    output:
        narrow_peak = "Results/07_peak/narrow_q{q_threshold}/{sample}_pe_peaks.narrowPeak"
    conda:
        config["macs_env"]
    group: "processing_group"
    params:
        gsize = config["deeptools"]["genome_sizes"][config["Spe"]],
        option = config["MACS"]["params"],
        outdir="Results/07_peak/narrow_q{q_threshold}/"
    log:
        log = "Results/07_peak/logs/{sample}_{q_threshold}.log"
    shell:
        """
        macs3 callpeak {params.option} -q {wildcards.q_threshold} -f BAMPE \
            -t {input.bam} -g {params.gsize} -n {wildcards.sample}_pe \
            --outdir {params.outdir} > {log.log} 2>&1
        """

rule callpeak_se:
    input:
        bam = get_uniq_list
    output:
        narrow_peak = "Results/07_peak/narrow_q{q_threshold}/{sample}_se_peaks.narrowPeak"
    conda:
        config["macs_env"]
    group: "processing_group"
    params:
        gsize = config["deeptools"]["genome_sizes"][config["Spe"]],
        option = config["MACS"]["params"],
        outdir="Results/07_peak/narrow_q{q_threshold}/"
    log:
        log = "Results/07_peak/logs/{sample}_{q_threshold}.log"
    shell:
        """
        macs3 callpeak {params.option} -q {wildcards.q_threshold} -f BAM \
            -t {input.bam} -g {params.gsize} -n {wildcards.sample}_se \
            --outdir {params.outdir} > {log.log} 2>&1
        """

rule extract_peak_pos:
    input: peak = "Results/07_peak/narrow_q{q_threshold}/{sample}_{dt}_peaks.narrowPeak"
    output:
        peak_pos = "Results/07_peak/narrow_q{q_threshold}/{sample}_{dt}_npks.bed"
    group: "processing_group"
    shell:
        """
        # 提取前三列（peak的坐标信息）
        cut -f 1-3 {input.peak} | sort -k1,1 -k2,2n > {output.peak_pos}
        """
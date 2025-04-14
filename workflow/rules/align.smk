rule align_reads_with_bowtie2_pe:
    input:
        read = get_trimmed_list
    output:
        sam = temp("Results/04_align/{sample}_pe.bowtie2_aln.sorted.sam")
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["bowtie2"]["params"],
        index= config["bowtie2"]["index"][config["Spe"]],
        align_out = "Results/04_align"
    log:
        log = "Results/04_align/logs/{sample}.log"
    shell:
        """
        bowtie2 {params.option} --no-mixed --no-discordant -x {params.index} -1 {input.read[0]} -2 {input.read[1]} -S {output.sam} > {log.log} 2>&1
        """

rule align_reads_with_bowtie2_se:
    input:
        read = get_trimmed_list
    output:
        sam = temp("Results/04_align/{sample}_se.bowtie2_aln.sorted.sam")
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["bowtie2"]["params"],
        index= config["bowtie2"]["index"][config["Spe"]],
        align_out = "Results/04_align"
    log:
        log = "Results/04_align/logs/{sample}.log"
    shell:
        """
        bowtie2 {params.option} -x {params.index} -U {input.read} -S {output.sam} > {log.log} 2>&1
        """

rule align_reads_with_bowtie2_human_pe:
    input:
        read = get_trimmed_list
    output:
        sam = temp("Results/04_align/{sample}_human_pe.bowtie2_aln.sorted.sam")
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["bowtie2"]["params"],
        index= config["bowtie2"]["index"]["human"],
        align_out = "Results/04_align"
    log:
        log = "Results/04_align/logs/{sample}_human.log"
    shell:
        """
        bowtie2 {params.option} --no-mixed --no-discordant -x {params.index} -1 {input.read[0]} -2 {input.read[1]} -S {output.sam} > {log.log} 2>&1
        """

rule align_reads_with_bowtie2_human_se:
    input:
        read = get_trimmed_list
    output:
        sam = temp("Results/04_align/{sample}_human_se.bowtie2_aln.sorted.sam")
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["bowtie2"]["params"],
        index= config["bowtie2"]["index"]["human"],
        align_out = "Results/04_align"
    log:
        log = "Results/04_align/logs/{sample}_human.log"
    shell:
        """
        bowtie2 {params.option} -x {params.index} -U {input.read} -S {output.sam} > {log.log} 2>&1
        """


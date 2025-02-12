rule align_reads_with_bowtie2:
    input:
        read = get_trimmed_list
    output:
        sam = temp("{align_out}/{sample}.bowtie2_aln.sorted.sam")
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        option = config["bowtie2"]["params"],
        index= config["bowtie2"]["index"][config["Spe"]],
        align_out = directories["align_out"]
    log:
        log = "{align_out}/logs/{sample}.log"
    shell:
        """
        if [ "{config[dt]}" == "SE" ]; then
           bowtie2 {params.option} -x {params.index} -U {input.read} -S {output.sam}
        else
           bowtie2 {params.option} --no-mixed --no-discordant -x {params.index} -1 {input.read[0]} -2 {input.read[1]} -S {output.sam}
        fi

        echo "alignment finished at $(date)"
        """


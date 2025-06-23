rule call_motif:
    input: narrow_peak = "Results/07_peak/narrow_q{q_threshold}/{sample}_{dt}_peaks.narrowPeak"
    output: "Results/07_peak/narrow_q{q_threshold}/motifs/{sample}_{dt}/knownResults.txt"
    conda:
        config["conda_env"]
    group: "Additional_analysis"
    params:
        genome=config["Spe"],
        option=config["find_motif"]["params"],
        output_dir= lambda wildcards: f"Results/07_peak/narrow_q{wildcards.q_threshold}/motifs/{wildcards.sample}_{wildcards.dt}"
    shell:
        """
        findMotifsGenome.pl {input.narrow_peak} {params.genome} {params.output_dir} {params.option} -p 20
        """
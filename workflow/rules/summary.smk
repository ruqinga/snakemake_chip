# visualization 结束之后运行summary
rule summary:
    input:
        callpeak_finish = [
            f"Results/07_peak/narrow_q{config['MACS']['q_threshold']}/{sample}_{'pe' if sample_info[sample] == 'PE' else 'se'}_npks.bed"
            for sample in samples
        ]
    output:
        summary = "Results/summary.csv"
    params:
        trim_log_dir="Results/02_trim_out/logs/",
        flagstat_dir="Results/05_dedu/stat/",
        align_logs_dir="Results/04_align/logs/",
        peak_dir = lambda wildcards: f"Results/07_peak/narrow_q{config['MACS']['q_threshold']}/"
    conda:
        config["conda_env"]
    group: "global_process"
    script:
        "../scripts/summary.py"
# visualization 结束之后运行summary
rule summary:
    input:
        extract_finish = expand("Results/06_visualization/{sample}_unique_{binsize}.bw", sample = samples, binsize = config["bamCoverage"]["binsize"])
    output:
        summary = "Results/summary.csv"
    params:
        trim_log_dir="Results/02_trim_out/logs/",
        flagstat_dir="Results/05_dedu/stat/",
        align_logs_dir="Results/04_align/logs/"
    conda:
        config["conda_env"]
    group: "global_process"
    script:
        "../scripts/summary.py"
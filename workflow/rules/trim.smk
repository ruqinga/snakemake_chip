import os

rule trim_pe:
    input:
        read = get_fq_list
    output:
        trimmed_read = [
                "Results/02_trim_out/{sample}_1_val_1.fq.gz",
                "Results/02_trim_out/{sample}_2_val_2.fq.gz"
        ]
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        rawdata = os.path.abspath(config["fq_dir"]),  # 将相对路径转换为绝对路径
        option = config["trim"]["params"],
        trim_out= "Results/02_trim_out"
    log:
        log = "Results/02_trim_out/logs/{sample}.log"
    shell:
        """
            # link input fq as Results/01_rawdata/
            ln -sfn {params.rawdata} Results/01_rawdata
            # trim_galore
            trim_galore {params.option} --paired {input.read[0]} {input.read[1]} -o {params.trim_out} > {log.log} 2>&1
        """

rule trim_se:
    input:
        read = get_fq_list
    output:
        trimmed_read = "Results/02_trim_out/{sample}_trimmed.fq.gz"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        rawdata = os.path.abspath(config["fq_dir"]),  # 将相对路径转换为绝对路径
        option = config["trim"]["params"],
        trim_out = "Results/02_trim_out"
    log:
        log = "Results/02_trim_out/logs/{sample}.log"
    shell:
        """
            # link input fq as Results/01_rawdata/
            ln -sfn {params.rawdata} Results/01_rawdata
            # trim_galore
            trim_galore {params.option} {input.read} -o {params.trim_out} > {log.log} 2>&1
        """
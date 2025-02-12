# 生成目录
def get_directories(config):
    """
    根据配置文件中的 work_dir 和目录结构返回所有目录的路径。
    """
    work_dir = config["work_dir"]
    directories = {
        key: value.format(work_dir=work_dir) if "{work_dir}" in value else value
        for key, value in config["directories"].items()
    }
    # 打印调试信息
    print(f"Resolved directories: {directories}")
    return directories

# 获取输入文件的basename（不带路径和扩展名）
def get_sample_name(filepath, dt):
    """
    从文件路径中提取样本名，基于数据类型（单端或双端）。
    """
    try:
        if dt == "PE":
            return filepath.split("/")[-1].replace("_1.fastq.gz", "").replace("_2.fastq.gz", "")
        elif dt == "SE":
            return filepath.split("/")[-1].replace(".fastq.gz", "")
        else:
            raise ValueError(f"Unsupported 'dt': {dt}")
    except Exception as e:
        raise ValueError(f"Error in get_sample_name for {filepath}: {e}")

def get_sample_list(config):
    """
    从配置文件中提取样本列表。
    """
    dt = config.get("dt")
    reads = config.get("reads", [])
    samples = {get_sample_name(read["read1"], dt) for read in reads}
    samples_list = list(samples)
    # 打印调试信息
    print("Samples:", samples_list)
    return samples_list

def get_all(directories, samples):
    """
    根据配置生成所有需要的目标文件路径列表。
    """
    dt = config.get("dt")
    is_pe = dt == "PE"

    binsize = config.get("bamCoverage", {}).get("binsize")

    # 定义不同模式下的文件扩展名
    trim_ext = "_1_val_1.fq.gz" if is_pe else "_trimmed.fq.gz"

    # 通用文件路径
    targets = [
        "{trim_out}/{sample}" + trim_ext,
        "{align_out}/{sample}.bowtie2_aln.sorted.sam",
        "{dedu_out}/{sample}_mapped_sorted_dedu_uniq.bam",
        "{dedu_out}/{sample}_mapped_sorted_dedu_uniq.bam.bai",
        "{bw_out}/{sample}_unique_{binsize}.bw",
        "{bw_out}/bed/{sample}_unique_{binsize}.bed"
    ]

    # 根据模式展开所有目标文件路径
    all_targets = [
        expand(path, **directories, sample=samples, binsize=binsize) for path in targets
    ]

    # 展平列表并打印调试信息
    flattened_targets = [item for sublist in all_targets for item in sublist]

    # 对bam文件进行stat和bai index
    suffixes = ["mapped_sorted", "mapped_sorted_dedu", "mapped_sorted_dedu_uniq"]
    flagstat_targets = expand(
        "{dedu_out}/stat/{sample}_{suffix}.flagstat",
        dedu_out=directories["dedu_out"], sample=samples, suffix=suffixes
    )

    # 合并所有目标路径
    all_files = flattened_targets + flagstat_targets

    print(all_files)  # 打印调试信息

    return all_files

def get_fq_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['fq_dir']}/{wildcards.sample}.fastq.gz"
    elif config["dt"] == "PE":
        return [
            f"{directories['fq_dir']}/{wildcards.sample}_1.fastq.gz",
            f"{directories['fq_dir']}/{wildcards.sample}_2.fastq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_trimmed_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['trim_out']}/{wildcards.sample}_trimmed.fq.gz"
    elif config["dt"] == "PE":
        return [
            f"{directories['trim_out']}/{wildcards.sample}_1_val_1.fq.gz",
            f"{directories['trim_out']}/{wildcards.sample}_2_val_2.fq.gz"
        ]
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_alined_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['align_out']}/{wildcards.sample}.bowtie2_aln.sorted.sam"
    elif config["dt"] == "PE":
        return f"{directories['align_out']}/{wildcards.sample}.bowtie2_aln.sorted.sam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

def get_uniq_list(wildcards):
    if config["dt"] == "SE":
        return f"{directories['dedu_out']}/{wildcards.sample}_mapped_sorted_dedu_uniq.bam"
    elif config["dt"] == "PE":
        return f"{directories['dedu_out']}/{wildcards.sample}_mapped_sorted_dedu_uniq.bam"
    else:
        raise ValueError(f"Invalid 'dt' configuration: {config['dt']}")

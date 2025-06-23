import re
import logging
import pandas as pd
from pathlib import Path
import sys

# trim_logs_dir = Path("Results/02_trim_out/logs/")
# flagstat_dir = Path("Results/05_dedu/stat/")
# align_logs_dir = Path("Results/04_align/logs/")
# peak_dir = Path("Results/07_peak/narrow_q0.05")
#
# outputfile = Path("Results/summary.txt")

# 获取 Snakemake 参数
trim_logs_dir = Path(snakemake.params.trim_log_dir)
flagstat_dir = Path(snakemake.params.flagstat_dir)
align_logs_dir = Path(snakemake.params.align_logs_dir)
peak_dir = Path(snakemake.params.peak_dir)

outputfile = Path(snakemake.output.summary)

# 设置日志格式
logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

def extract_from_trim_log(log_file):
    with open(log_file, 'r') as f:
        content = f.read()

    # Primary regex patterns
    input_reads = re.search(r"Total number of sequences analysed:\s*(\d+)", content)
    less_than_cutoff = re.search(r"Number of sequence pairs removed because at least one read was shorter than the length cutoff \(30 bp\):\s*(\d+)", content)
    N_than_4 = re.search(r"Number of sequence pairs removed because at least one read contained more N\(s\) than the specified limit of 4:\s*(\d+)", content)

    # Secondary fallback patterns
    if not input_reads or not less_than_cutoff or not N_than_4:
        input_reads = re.search(r"(\d+) sequences processed in total", content)
        less_than_cutoff = re.search(r"Sequences removed because they became shorter than the length cutoff of 30 bp:\s*(\d+)", content)
        N_than_4 = re.search(r"Sequences removed because they contained more Ns than the cutoff of 4:\s*(\d+)", content)

    # Convert matches to integers (default to 0 if not found)
    input_reads_value = int(input_reads.group(1)) if input_reads else 0
    less_than_cutoff_value = int(less_than_cutoff.group(1)) if less_than_cutoff else 0
    N_than_4_value = int(N_than_4.group(1)) if N_than_4 else 0

    # Logging warnings if values are missing
    if not input_reads:
        logging.warning(f"{log_file} 没有找到 input_reads")
    if not less_than_cutoff:
        logging.warning(f"{log_file} 没有找到 less_than_cutoff")
    if not N_than_4:
        logging.warning(f"{log_file} 没有找到 N_than_4")

    # Calculate trimmed reads
    trimmed_reads = input_reads_value - less_than_cutoff_value - N_than_4_value
    return input_reads_value, trimmed_reads

def extract_from_align_log(log_file):
    """Extract aligned read counts from a bowtie2 log file."""
    with open(log_file, 'r') as f:
        content = f.read()

    # Primary pattern: "aligned concordantly exactly 1 time" and ">1 times"
    exactly_1_match = re.search(r"(\d+) \(.*?\) aligned concordantly exactly 1 time", content)
    more_than_1_match = re.search(r"(\d+) \(.*?\) aligned concordantly >1 times", content)

    # Secondary pattern: "aligned exactly 1 time" and ">1 times"
    if not exactly_1_match or not more_than_1_match:
        exactly_1_match = re.search(r"(\d+) \(.*?\) aligned exactly 1 time", content)
        more_than_1_match = re.search(r"(\d+) \(.*?\) aligned >1 times", content)

    # Convert matches to integers and sum them
    aligned_reads = sum(int(m.group(1)) for m in [exactly_1_match, more_than_1_match] if m)
    return aligned_reads

def extract_from_flagstat(flagstat):
    with open(flagstat, 'r') as f:
        content = f.read()

    reads = re.search(r"(\d+)\s+\+\s+\d+\s+in total (QC-passed reads + QC-failed reads)", content)
    #  samtools view -c CT_3CHA09-HA_mapped_sorted_dedu_uniq.bam

    if not reads:
        logging.error(f"{flagstat} 缺失配对相关字段")
        sys.exit(1)

    return int(reads.group(1))

def count_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    if len(lines) == 0:
        logging.warning(f"{file_path} 没有找到peak")

    return len(lines)

def main():
    raw_reads_dict = {}
    clean_reads_dict = {}
    mapped_reads_dict = {}
    mapping_rate_dict = {}
    mapping_rate_dict_human = {}
    dedu_reads_dict = {}
    duplicates_rate_dict = {}
    uniq_reads_dict = {}
    uniq_ratio_dict = {}
    peaknum_dic = {}


    # 1. trim log
    for log_file in trim_logs_dir.glob("*.log"):
        sample = log_file.stem
        raw, clean = extract_from_trim_log(log_file)
        raw_reads_dict[sample] = raw
        clean_reads_dict[sample] = clean

    # 2. flagstat
    for log_file in flagstat_dir.glob("*_mapped_sorted.flagstat"):
        sample = re.sub(r"_mapped_sorted\.flagstat$", "", log_file.name)
        mapped = extract_from_flagstat(log_file)
        mapped_reads_dict[sample] = mapped
        print(f"Sample: {sample}, Mapped Reads: {mapped}")
        mapping_rate_dict[sample] = mapped / clean_reads_dict.get(sample, 1)

    for log_file in flagstat_dir.glob("*_mapped_sorted_dedu.flagstat"):
        sample = re.sub(r"_mapped_sorted_dedu\.flagstat$", "", log_file.name)
        dedu = extract_from_flagstat(log_file)
        dedu_reads_dict[sample] = dedu
        print(f"Sample: {sample}, dedu Reads: {mapped}")
        duplicates_rate_dict[sample] =(mapped_reads_dict.get(sample, 1) - dedu) / mapped_reads_dict.get(sample, 1)

    for log_file in flagstat_dir.glob("*_mapped_sorted_dedu_uniq.flagstat"):
        sample = re.sub(r"_mapped_sorted_dedu_uniq\.flagstat$", "", log_file.name)
        uniq = extract_from_flagstat(log_file)
        uniq_reads_dict[sample] = uniq
        uniq_ratio_dict[sample] = uniq / dedu_reads_dict.get(sample, 1)

    # 3. align log of human
    for log_file in align_logs_dir.glob("*_human*.log"):
        sample = re.sub(r"_human(_se|_pe)\.log$", "", log_file.name)
        #print(log_file.name)
        aligned_reads = extract_from_align_log(log_file)
        mapping_rate_dict_human[sample] = aligned_reads / clean_reads_dict.get(sample, 1)

    # 4. peak number
    for peakfile in peak_dir.glob("*_npks.bed"):
        sample = re.sub(r"(_se|_pe)_npks$", "", peakfile.stem)
        peak_number = count_lines(peakfile)

        peaknum_dic[sample] = peak_number

    # 6. 构建 DataFrame
    df = pd.DataFrame({
        "raw_reads": raw_reads_dict,
        "clean_reads": clean_reads_dict,
        "mapped_reads": mapped_reads_dict,
        "mapped_rate": mapping_rate_dict,
        "mapped_rate_human": mapping_rate_dict_human,
        "dedu_reads": dedu_reads_dict,
        "duplicates_rate": duplicates_rate_dict,
        "uniq_reads": uniq_reads_dict,
        "uniq_ratio": uniq_ratio_dict,
        "peaknum": peaknum_dic
    }).reset_index().rename(columns={"index": "filename"})

    # 格式化：整数列与浮点列分别处理
    int_cols = ["raw_reads", "clean_reads", "mapped_reads", "dedu_reads", "uniq_reads", "peaknum"]
    float_cols = ["mapped_rate", "mapped_rate_human", "duplicates_rate", "uniq_ratio"]

    for col in int_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0).astype(int)

    for col in float_cols:
        if col in df.columns:
            df[col] = df[col].astype(float).round(3)

    # 输出结果
    df.to_csv(outputfile, sep='\t', index=False, encoding='utf-8')
    print("Wrote summary file to {}".format(outputfile))


if __name__ == "__main__":
    main()

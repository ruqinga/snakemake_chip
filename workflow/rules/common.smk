# 解析input，生成output
class SampleProcessor:
    def __init__(self, config):
        """
        初始化 SampleProcessor。

        参数:
            config (dict): 包含 reads 和 bamCoverage 的配置信息。
        """
        self.reads = config.get("reads", [])
        self.sample_names = []
        self.sample_info = {}
        self.binsize = config.get("bamCoverage", {}).get("binsize", 100)
        self.q_threshold = config.get("MACS", {}).get("q_threshold", 0.05)
        self.process_reads()

    def process_reads(self):
        """
        处理 reads 列表，提取样本名并确定是单端 (SE) 还是双端 (PE)。
        """
        for read in self.reads:
            sample_name = read["read1"].split("/")[-1].replace("_1.fastq.gz", "").replace(".fastq.gz", "")
            if sample_name in self.sample_names:
                raise ValueError(f"Duplicate sample name '{sample_name}' found.")
            self.sample_names.append(sample_name)

            if "read2" in read and read["read2"]:
                self.sample_info[sample_name] = "PE"
            else:
                self.sample_info[sample_name] = "SE"

    def get_trim_ext(self, sample):
        """
        根据样本类型返回 trim 后的文件后缀名。

        参数:
            sample (str): 样本名。

        返回:
            str: PE 返回 "1_val_1"，SE 返回 "trimmed"
        """
        return "1_val_1" if self.sample_info[sample] == "PE" else "trimmed"

    def generate_targets(self, sample, binsize,q_threshold):
        """
        生成给定样本的所有目标文件路径。

        参数:
            sample (str): 样本名
            binsize (int): 可视化使用的 bin size

        返回:
            list[str]: 目标路径列表
        """
        dt = 'pe' if self.sample_info[sample] == 'PE' else 'se'
        suffixes = ["mapped_sorted", "mapped_sorted_dedu", "mapped_sorted_dedu_uniq"]
        return [
            f"Results/02_trim_out/{sample}_{self.get_trim_ext(sample)}.fq.gz",
            f"Results/03_qc/rawdata/multiqc_report.html",
            f"Results/03_qc/cleandata/multiqc_report.html",
            f"Results/04_align/logs/{sample}_human.log",
            f"Results/05_dedu/{sample}_mapped_sorted_dedu_uniq.bam",
            *[f"Results/05_dedu/stat/{sample}_{s}.flagstat" for s in suffixes],
            f"Results/06_visualization/{sample}_unique_{binsize}.bw",
            f"Results/06_visualization/bed/{sample}_unique_{binsize}.bed",
            f"Results/07_peak/narrow_q{q_threshold}/{sample}_{dt}_peaks.narrowPeak",
            f"Results/summary.csv"
        ]

    def get_all_targets(self):
        """
        获取所有样本的目标文件路径。

        返回:
            list[str]: 所有样本的目标路径列表
        """
        all_paths = []
        for sample in self.sample_names:
            paths = self.generate_targets(sample, self.binsize,self.q_threshold)
            all_paths.extend(paths)
        return all_paths



def get_fq_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"{config['fq_dir']}/{wildcards.sample}.fastq.gz"
    else:
        return [
            f"{config['fq_dir']}/{wildcards.sample}_1.fastq.gz",
            f"{config['fq_dir']}/{wildcards.sample}_2.fastq.gz"
        ]


def get_trimmed_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/02_trim_out/{wildcards.sample}_trimmed.fq.gz"
    else:
        return [
            f"Results/02_trim_out/{wildcards.sample}_1_val_1.fq.gz",
            f"Results/02_trim_out/{wildcards.sample}_2_val_2.fq.gz"
        ]

def get_alined_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/04_align/{wildcards.sample}_se.bowtie2_aln.sorted.sam"
    else:
        return f"Results/04_align/{wildcards.sample}_pe.bowtie2_aln.sorted.sam"


def get_uniq_list(wildcards):
    if sample_info[wildcards.sample] == "SE":
        return f"Results/05_dedu/{wildcards.sample}_mapped_sorted_dedu_uniq.bam"
    else:
        return f"Results/05_dedu/{wildcards.sample}_mapped_sorted_dedu_uniq.bam"


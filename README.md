# snakemake_chip
todo list:
- 对log文件夹加上保护，避免在删除整个文件夹时被删除

用法: bash run.sh <fq_dir> -y
- <fq_dir> 为fastq.gz文件的存放目录，单双端以 _1.fastq.gz _2.fastq.gz 区分
- -y 表示实际运行，否则只会打印出要运行的命令

单独提交某个任务
```shell
fq_dir="../rawdata"
reads_json=$(cat "$fq_dir/reads_json.json")

snakemake \
    -np \
    --use-conda \
    --config fq_dir="$fq_dir" reads="$reads_json" \
    --allowed-rules call_motif

# --touch --ignore-incomplete --cleanup-metadata 都是方法之一
snakemake --touch -np --use-conda --config fq_dir="$fq_dir" reads="$reads_json" -f --until call_motif --ignore-incomplete
snakemake --cleanup-metadata Results/07_peak/narrow_q0.05/*.narrowPeak
```
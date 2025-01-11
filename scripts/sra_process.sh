# sra处理函数
sra_process() {
    local download=$1
    local srrid=$2
    local expected_md5=$3
    local rename=$4
    local sra_file="${sra_dir}/${srrid}/${srrid}.sra"

    # 切换到目标目录
    cd "$sra_dir" || { echo "无法切换到目录 $sra_dir"; return 1; }
    
    # 使用 prefetch 下载 SRA 文件
    if download=True
    echo "downloading ${srrid}"
    prefetch --max-size 50G "$srrid" || { echo "prefetch 下载失败"; return 1; }
    else
      echo "$sra_file already exist"

    
    # MD5 校验
    calculated_md5=$(md5sum "$sra_file" | awk '{print $1}')
    if [[ $calculated_md5 != $expected_md5 ]]; then
        echo "MD5 校验和不匹配，下载可能出错。"
        return 1
    fi
    echo "MD5 校验通过。"

    # 重命名 SRA 文件
    mv "${sra_dir}/${srrid}/${srrid}.sra" "${sra_dir}/${srrid}/${rename}.sra" || { echo "重命名失败"; return 1; }

    # 解压为 FASTQ 格式
    fasterq-dump --threads 10 --split-3 --outfile "${fq_dir}/${rename}.fastq" "${sra_dir}/${srrid}/${rename}.sra" || { echo "文件 "${sra_dir}/${srrid}/${rename}.sra" 解压失败。"; return 1; }

    echo "文件 "${sra_dir}/${srrid}/${rename}.sra" 解压成功。"
    }

expert -f sra_process

sra_dir = ""
export sra_dir
parallel --header : --colsep '\t' --link sra_process "True" {accession} {sra_md5} {sample_title}  :::: name_srr.txt

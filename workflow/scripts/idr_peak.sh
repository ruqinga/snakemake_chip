#!/bin/bash
#nohup bash idr.sh > idr.log 2>&1 &

conda activate idr
nPeak_dir="/home_data/home/slst/leixy2023/data/project/CTCF/20250519_peak_collect/peak"
output_dir="/home_data/home/slst/leixy2023/data/project/CTCF/20250519_peak_collect/idr"

# 提取同一sample的不同rep
## 写for循环而不是while的原因：while 循环内部的变量赋值不会影响外部的变量，除非改的麻烦点。
## 创建关联数组
declare -A sample_reps
## 清除数组 unset sample_reps
## 循环
files=$(ls ${nPeak_dir}/*.narrowPeak)
for file in $files; do
    # 提取样本名（假设样本名是文件名中 "rep" 之前的部分）
    basename=$(basename "$file")
    sample=$(echo "$basename" | sed -E 's/_rep[0-9]+_.*//')
    # 将文件名添加到对应的样本数组中
    sample_reps[$sample]+="$file "
done
# 检测 echo ${sample_reps[$sample]}


# idr
for sample in "${!sample_reps[@]}"; do
    idr_output="${output_dir}/${sample}_idr.narrowPeak"

    # 执行 IDR 分析
    idr --samples ${sample_reps[$sample]} \
        --input-file-type narrowPeak \
        --rank p.value \
        --output-file "$idr_output" \
        --plot
    echo "✓ $sample IDR finished for: ${sample_reps[$sample]} -> $idr_output"
done

## 后续如果要合并超过两个rep的话，可以通过下面的方法提取
#IFS=' ' read -r -a files <<< "${sample_reps[$sample]}"
## 确保有且仅有两个文件路径
#if [ "${#files[@]}" -ne 2 ]; then
#    echo "Error: Sample $sample does not have exactly two replicates."
#    continue
#fi
#
#idr --samples "${files[0]}" "${files[1]}"


rule dedu_and_uniq:
    input:
        sam = get_alined_list
    output:
        mapped_bam = temp("{dedu_out}/{sample}_mapped.bam"),
        sorted_bam = "{dedu_out}/{sample}_mapped_sorted.bam",
        dedu_bam = "{dedu_out}/{sample}_mapped_sorted_dedu.bam",
        uniq_bam = "{dedu_out}/{sample}_mapped_sorted_dedu_uniq.bam",
        uniq_bam_bai = "{dedu_out}/{sample}_mapped_sorted_dedu_uniq.bam.bai"
    conda:
        config["conda_env"]
    group: "processing_group"
    params:
        dedu_out = directories["dedu_out"],
        thread = config["threads"]
    log:
        log = "{dedu_out}/logs/{sample}.log"
    shell:
        """
        # 去掉未比对的序列和其配对序列
        sambamba view -t {params.thread} -F "not (unmapped or mate_is_unmapped)" -f bam -S {input.sam} -o {output.mapped_bam} > {log.log} 2>&1
        sambamba sort {output.mapped_bam} -o {output.sorted_bam} > {log.log} 2>&1
        
        # 去除 PCR 扩增产生的重复序列
        sambamba markdup -t {params.thread} -r {output.sorted_bam} {output.dedu_bam} > {log.log} 2>&1
        
        # 通过 grep -v "XS:i:" 过滤掉带次优比对标记的 reads，只保留唯一比对（unique alignment）的序列。
        samtools view -Sh -@ {params.thread} {output.dedu_bam} | grep -v "XS:i:" | samtools sort -O bam -@ {params.thread} -o {output.uniq_bam} > {log.log} 2>&1
        sambamba index {output.uniq_bam} > {output.uniq_bam_bai}
        
        echo "dedu finished at $(date)"
        """

rule samtools_flagstat:
    input:
        bam="{dedu_out}/{sample}_{suffix}.bam"
    output:
        flagstat="{dedu_out}/stat/{sample}_{suffix}.flagstat"
    conda:
        config["conda_env"]
    group: "processing_group"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """

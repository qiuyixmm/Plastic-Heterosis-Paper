#!/usr/bin/env bash

# 输出文件
out="hybrid_OV_matched_summary.tsv"

# 表头
echo -e "sample\tPeking_reads\tPeking_matched_reads\tPeking_matched_ratio\tMuscovy_reads\tMuscovy_matched_reads\tMuscovy_matched_ratio" > "$out"

# 循环所有 hybrid featureCounts 文件
for hybrid in \
Pek.OV.13.Mus.OV.16.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.17.Mus.OV.20.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.1.Mus.OV.32.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.21.Mus.OV.4.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.25.Mus.OV.36.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.29.Mus.OV.8.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.33.Mus.OV.40.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.37.Mus.OV.28.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.5.Mus.OV.24.hisat2.filter.sorted.bam.featureCounts \
Pek.OV.9.Mus.OV.12.hisat2.filter.sorted.bam.featureCounts
do
    # sample 名（去掉后缀）
    sample=${hybrid%.hisat2.filter.sorted.bam.featureCounts}

    # 拆出 Pek 和 Mus 的单亲文件名
    Pek_file=$(echo "$sample" | sed -E 's/(Pek\.OV\.[0-9]+).*/\1.hisat2.filter.sorted.bam.featureCounts/')
    Mus_file=$(echo "$sample" | sed -E 's/.*(Mus\.OV\.[0-9]+)/\1.hisat2.filter.sorted.bam.featureCounts/')

    # overlap 文件
    Pek_overlap="${sample}.Pek.overlap.featureCounts"
    Mus_overlap="${sample}.Mus.overlap.featureCounts"

    # 1. 计算 overlap（整行一致）
    awk 'NR==FNR{a[$0]; next} $0 in a' "$Pek_file" "$hybrid" > "$Pek_overlap"
    awk 'NR==FNR{a[$0]; next} $0 in a' "$Mus_file" "$hybrid" > "$Mus_overlap"

    # 2. 统计 reads
    Peking_reads=$(awk '$1~/^Pek/ && $2=="Assigned"{c++} END{print c+0}' "$hybrid")
    Peking_matched_reads=$(awk '$1~/^Pek/ && $2=="Assigned"{c++} END{print c+0}' "$Pek_overlap")

    Muscovy_reads=$(awk '$1~/^Mus/ && $2=="Assigned"{c++} END{print c+0}' "$hybrid")
    Muscovy_matched_reads=$(awk '$1~/^Mus/ && $2=="Assigned"{c++} END{print c+0}' "$Mus_overlap")

    # 3. 比例（浮点）
    Peking_ratio=$(awk -v m="$Peking_matched_reads" -v t="$Peking_reads" \
        'BEGIN{if(t>0) printf "%.6f", m/t; else print "NA"}')

    Muscovy_ratio=$(awk -v m="$Muscovy_matched_reads" -v t="$Muscovy_reads" \
        'BEGIN{if(t>0) printf "%.6f", m/t; else print "NA"}')

    # 4. 写入结果
    echo -e "${sample}\t${Peking_reads}\t${Peking_matched_reads}\t${Peking_ratio}\t${Muscovy_reads}\t${Muscovy_matched_reads}\t${Muscovy_ratio}" >> "$out"
done
make_summary.sh (END)

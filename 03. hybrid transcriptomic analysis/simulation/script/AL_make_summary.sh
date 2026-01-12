#!/usr/bin/env bash

# 输出文件
out="hybrid_AL_matched_summary.tsv"

# 表头
echo -e "sample\tPeking_reads\tPeking_matched_reads\tPeking_matched_ratio\tMuscovy_reads\tMuscovy_matched_reads\tMuscovy_matched_ratio" > "$out"

# 循环所有 hybrid featureCounts 文件
for hybrid in \
Pek.AL.101.Mus.AL.96.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.105.Mus.AL.88.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.109.Mus.AL.112.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.113.Mus.AL.116.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.117.Mus.AL.120.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.81.Mus.AL.108.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.85.Mus.AL.84.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.93.Mus.AL.92.hisat2.filter.sorted.bam.featureCounts \
Pek.AL.97.Mus.AL.104.hisat2.filter.sorted.bam.featureCounts
do
    # sample 名（去掉后缀）
    sample=${hybrid%.hisat2.filter.sorted.bam.featureCounts}

    # 拆出 Pek 和 Mus 的单亲文件名
    Pek_file=$(echo "$sample" | sed -E 's/(Pek\.AL\.[0-9]+).*/\1.hisat2.filter.sorted.bam.featureCounts/')
    Mus_file=$(echo "$sample" | sed -E 's/.*(Mus\.AL\.[0-9]+)/\1.hisat2.filter.sorted.bam.featureCounts/')

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

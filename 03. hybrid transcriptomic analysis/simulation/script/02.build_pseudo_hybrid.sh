#!/bin/bash
# ====================================================================================================================
# 功能：从北京鸭和番鸭RNA-seq数据中随机组合样本，各抽取一半pair-end reads（保持配对），生成混合样本fastq.gz文件。
# 依赖：seqtk、seqkit、gzip
# ====================================================================================================================

set -euo pipefail

# -------------------------------
# 检查参数
# -------------------------------
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <PEK_DIR> <MUS_DIR> <OUT_DIR>"
    echo "Example: $0 /path/to/Pek.AL /path/to/Mus.AL /path/to/output"
    exit 1
fi

PEK_DIR="$1"
MUS_DIR="$2"
OUT_DIR="$3"
mkdir -p "${OUT_DIR}"

# 指定 seqkit 路径
SEQKIT="/gpfs/home/xumingmin/dir.xmm/bin/seqkit/seqkit"

# ----------------------------------------------------------
# 自动获取样本编号
# ----------------------------------------------------------
#PEK_SAMPLES=($(ls ${PEK_DIR}/Pek.AL.*_1.fastq.gz | sed 's/.*\/\(Pek\.AL\.[0-9]*\)_1\.fastq\.gz/\1/' | sort))
#MUS_SAMPLES=($(ls ${MUS_DIR}/Mus.AL.*_1.fastq.gz | sed 's/.*\/\(Mus\.AL\.[0-9]*\)_1\.fastq\.gz/\1/' | sort))

PEK_SAMPLES=($(ls ${PEK_DIR}/Pek.*_1.fastq.gz | xargs -n1 basename | sed 's/_1\.fastq\.gz//' | sort))
MUS_SAMPLES=($(ls ${MUS_DIR}/Mus.*_1.fastq.gz | xargs -n1 basename | sed 's/_1\.fastq\.gz//' | sort))

NUM_PEK=${#PEK_SAMPLES[@]}
NUM_MUS=${#MUS_SAMPLES[@]}
NUM_PAIRS=$(( NUM_PEK < NUM_MUS ? NUM_PEK : NUM_MUS ))

echo " 北京鸭样本数：${NUM_PEK}"
echo " 番鸭样本数：${NUM_MUS}"
echo " 将生成混合样本数量：${NUM_PAIRS}"

# 打乱样本顺序
shuf_peks=($(printf "%s\n" "${PEK_SAMPLES[@]}" | shuf))
shuf_muss=($(printf "%s\n" "${MUS_SAMPLES[@]}" | shuf))

# 截取最少一方的样本数
shuf_peks=("${shuf_peks[@]:0:$NUM_PAIRS}")
shuf_muss=("${shuf_muss[@]:0:$NUM_PAIRS}")

# ----------------------------------------------------------
# 开始生成混合样本
# ----------------------------------------------------------
echo "开始生成混合样本..."
echo "---------------------------------------"

for i in $(seq 0 $((NUM_PAIRS-1))); do
    pek=${shuf_peks[$i]}
    mus=${shuf_muss[$i]}

    pek_1="${PEK_DIR}/${pek}_1.fastq.gz"
    pek_2="${PEK_DIR}/${pek}_2.fastq.gz"
    mus_1="${MUS_DIR}/${mus}_1.fastq.gz"
    mus_2="${MUS_DIR}/${mus}_2.fastq.gz"

    out_prefix="${OUT_DIR}/${pek}.${mus}"
    echo "正在处理配对：${pek} + ${mus}"

    # ---------------------------
    # 北京鸭：抽取一半 reads 并保持配对
    # ---------------------------
    $SEQKIT sample -s 123 "${pek_1}" -p 0.5 | $SEQKIT seq --name --only-id > "${OUT_DIR}/pek_ids.txt"
    zcat "${pek_1}" | $SEQKIT grep -f "${OUT_DIR}/pek_ids.txt" > "${OUT_DIR}/pek_tmp_1.fastq"
    zcat "${pek_2}" | $SEQKIT grep -f "${OUT_DIR}/pek_ids.txt" > "${OUT_DIR}/pek_tmp_2.fastq"

    # ---------------------------
    # 番鸭：抽取一半 reads 并保持配对
    # ---------------------------
    $SEQKIT sample -s 456 "${mus_1}" -p 0.5 | $SEQKIT seq --name --only-id > "${OUT_DIR}/mus_ids.txt"
    zcat "${mus_1}" | $SEQKIT grep -f "${OUT_DIR}/mus_ids.txt" > "${OUT_DIR}/mus_tmp_1.fastq"
    zcat "${mus_2}" | $SEQKIT grep -f "${OUT_DIR}/mus_ids.txt" > "${OUT_DIR}/mus_tmp_2.fastq"

    # ---------------------------
    # 合并北京鸭 + 番鸭 reads 并压缩
    # ---------------------------
    cat "${OUT_DIR}/pek_tmp_1.fastq" "${OUT_DIR}/mus_tmp_1.fastq" | gzip > "${out_prefix}_1.fastq.gz"
    cat "${OUT_DIR}/pek_tmp_2.fastq" "${OUT_DIR}/mus_tmp_2.fastq" | gzip > "${out_prefix}_2.fastq.gz"

    # 清理临时文件
    rm "${OUT_DIR}/pek_tmp_"*.fastq "${OUT_DIR}/mus_tmp_"*.fastq "${OUT_DIR}/pek_ids.txt" "${OUT_DIR}/mus_ids.txt"
done

echo "---------------------------------------"
echo "混合样本生成完成！输出目录：${OUT_DIR}"


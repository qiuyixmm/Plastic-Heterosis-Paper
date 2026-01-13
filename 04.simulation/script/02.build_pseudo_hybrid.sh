#!/bin/bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <PEK_DIR> <MUS_DIR> <OUT_DIR>"
    echo "Example: $0 /path/to/Pek.AL /path/to/Mus.AL /path/to/output"
    exit 1
fi

PEK_DIR="$1"
MUS_DIR="$2"
OUT_DIR="$3"
mkdir -p "${OUT_DIR}"

SEQKIT="/gpfs/home/xumingmin/dir.xmm/bin/seqkit/seqkit"

# ----------------------------------------------------------
# obtain sample ID
# ----------------------------------------------------------
#PEK_SAMPLES=($(ls ${PEK_DIR}/Pek.AL.*_1.fastq.gz | sed 's/.*\/\(Pek\.AL\.[0-9]*\)_1\.fastq\.gz/\1/' | sort))
#MUS_SAMPLES=($(ls ${MUS_DIR}/Mus.AL.*_1.fastq.gz | sed 's/.*\/\(Mus\.AL\.[0-9]*\)_1\.fastq\.gz/\1/' | sort))

PEK_SAMPLES=($(ls ${PEK_DIR}/Pek.*_1.fastq.gz | xargs -n1 basename | sed 's/_1\.fastq\.gz//' | sort))
MUS_SAMPLES=($(ls ${MUS_DIR}/Mus.*_1.fastq.gz | xargs -n1 basename | sed 's/_1\.fastq\.gz//' | sort))

NUM_PEK=${#PEK_SAMPLES[@]}
NUM_MUS=${#MUS_SAMPLES[@]}
NUM_PAIRS=$(( NUM_PEK < NUM_MUS ? NUM_PEK : NUM_MUS ))

echo " Peking ducks：${NUM_PEK}"
echo " Nuscovy ducks：${NUM_MUS}"
echo " Hybrid ducks：${NUM_PAIRS}"

# resample
shuf_peks=($(printf "%s\n" "${PEK_SAMPLES[@]}" | shuf))
shuf_muss=($(printf "%s\n" "${MUS_SAMPLES[@]}" | shuf))

shuf_peks=("${shuf_peks[@]:0:$NUM_PAIRS}")
shuf_muss=("${shuf_muss[@]:0:$NUM_PAIRS}")

for i in $(seq 0 $((NUM_PAIRS-1))); do
    pek=${shuf_peks[$i]}
    mus=${shuf_muss[$i]}

    pek_1="${PEK_DIR}/${pek}_1.fastq.gz"
    pek_2="${PEK_DIR}/${pek}_2.fastq.gz"
    mus_1="${MUS_DIR}/${mus}_1.fastq.gz"
    mus_2="${MUS_DIR}/${mus}_2.fastq.gz"

    out_prefix="${OUT_DIR}/${pek}.${mus}"

    # ---------------------------
    # Peking ducks: 50% downsampling
    # ---------------------------
    $SEQKIT sample -s 123 "${pek_1}" -p 0.5 | $SEQKIT seq --name --only-id > "${OUT_DIR}/pek_ids.txt"
    zcat "${pek_1}" | $SEQKIT grep -f "${OUT_DIR}/pek_ids.txt" > "${OUT_DIR}/pek_tmp_1.fastq"
    zcat "${pek_2}" | $SEQKIT grep -f "${OUT_DIR}/pek_ids.txt" > "${OUT_DIR}/pek_tmp_2.fastq"

    # ---------------------------
    # Muscovy ducks: 50% downsampling
    # ---------------------------
    $SEQKIT sample -s 456 "${mus_1}" -p 0.5 | $SEQKIT seq --name --only-id > "${OUT_DIR}/mus_ids.txt"
    zcat "${mus_1}" | $SEQKIT grep -f "${OUT_DIR}/mus_ids.txt" > "${OUT_DIR}/mus_tmp_1.fastq"
    zcat "${mus_2}" | $SEQKIT grep -f "${OUT_DIR}/mus_ids.txt" > "${OUT_DIR}/mus_tmp_2.fastq"

    # ---------------------------
    # merge
    # ---------------------------
    cat "${OUT_DIR}/pek_tmp_1.fastq" "${OUT_DIR}/mus_tmp_1.fastq" | gzip > "${out_prefix}_1.fastq.gz"
    cat "${OUT_DIR}/pek_tmp_2.fastq" "${OUT_DIR}/mus_tmp_2.fastq" | gzip > "${out_prefix}_2.fastq.gz"


    rm "${OUT_DIR}/pek_tmp_"*.fastq "${OUT_DIR}/mus_tmp_"*.fastq "${OUT_DIR}/pek_ids.txt" "${OUT_DIR}/mus_ids.txt"
done

echo "---------------------------------------"
echo "done：${OUT_DIR}"


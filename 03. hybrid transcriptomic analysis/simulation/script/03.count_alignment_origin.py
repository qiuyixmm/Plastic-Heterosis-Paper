#!/gpfs/home/xumingmin/dir.xmm/bin/python3/bin/python
# ==================================================================================================
# 功能：统计混合样本（骡鸭）比对到 hybrid genome 的来源 reads 数量，并输出跨物种 BAM 文件
# 北京鸭染色体前缀：NC_ 或 NW_
# 番鸭染色体前缀：HiC_
# 输出：
#   1. 四类 reads 统计（Pek_to_Pek, Pek_to_Mus, Mus_to_Pek, Mus_to_Mus）
#   2. Pek_to_Mus.bam
#   3. Mus_to_Pek.bam
# ==================================================================================================

import pysam
import sys
import os
import subprocess

if len(sys.argv) != 2:
    print("Usage: count_alignment_origin.py <sorted.bam>")
    sys.exit(1)

bam_path = sys.argv[1]
sample = os.path.basename(bam_path).replace(".sorted.bam", "")

# 检查索引文件是否存在
if not os.path.exists(bam_path + ".bai"):
    print(f"[Info] BAM index not found for {bam_path}, creating index...", file=sys.stderr)
    subprocess.run(["samtools", "index", bam_path], check=True)

bam = pysam.AlignmentFile(bam_path, "rb")

# 输出跨物种 BAM 文件
pek_to_mus_bam = f"{sample}_Pek_to_Mus.bam"
mus_to_pek_bam = f"{sample}_Mus_to_Pek.bam"
pek_to_mus_out = pysam.AlignmentFile(pek_to_mus_bam, "wb", template=bam)
mus_to_pek_out = pysam.AlignmentFile(mus_to_pek_bam, "wb", template=bam)

# 初始化计数
counts = {
    "Pek_to_Pek": 0,
    "Pek_to_Mus": 0,
    "Mus_to_Pek": 0,
    "Mus_to_Mus": 0
}

# 遍历 BAM
for read in bam.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    qname = read.query_name
    rname = bam.get_reference_name(read.reference_id)

    # Pek reads
    if qname.startswith("Pek"):
        if rname.startswith("NC_") or rname.startswith("NW_"):
            counts["Pek_to_Pek"] += 1
        elif rname.startswith("HiC_"):
            counts["Pek_to_Mus"] += 1
            pek_to_mus_out.write(read)

    # Mus reads
    elif qname.startswith("Mus"):
        if rname.startswith("NC_") or rname.startswith("NW_"):
            counts["Mus_to_Pek"] += 1
            mus_to_pek_out.write(read)
        elif rname.startswith("HiC_"):
            counts["Mus_to_Mus"] += 1

bam.close()
pek_to_mus_out.close()
mus_to_pek_out.close()

# 为跨物种 BAM 建索引
subprocess.run(["samtools", "index", pek_to_mus_bam], check=True)
subprocess.run(["samtools", "index", mus_to_pek_bam], check=True)

# 输出统计结果
print(f"{sample}\t{counts['Pek_to_Pek']}\t{counts['Pek_to_Mus']}\t{counts['Mus_to_Pek']}\t{counts['Mus_to_Mus']}")

#!/gpfs/home/xumingmin/dir.xmm/bin/python3/bin/python
import pysam
import sys
import os
import subprocess

if len(sys.argv) != 2:
    print("Usage: count_alignment_origin.py <sorted.bam>")
    sys.exit(1)

bam_path = sys.argv[1]
sample = os.path.basename(bam_path).replace(".sorted.bam", "")

if not os.path.exists(bam_path + ".bai"):
    print(f"[Info] BAM index not found for {bam_path}, creating index...", file=sys.stderr)
    subprocess.run(["samtools", "index", bam_path], check=True)

bam = pysam.AlignmentFile(bam_path, "rb")

pek_to_mus_bam = f"{sample}_Pek_to_Mus.bam"
mus_to_pek_bam = f"{sample}_Mus_to_Pek.bam"
pek_to_mus_out = pysam.AlignmentFile(pek_to_mus_bam, "wb", template=bam)
mus_to_pek_out = pysam.AlignmentFile(mus_to_pek_bam, "wb", template=bam)

counts = {
    "Pek_to_Pek": 0,
    "Pek_to_Mus": 0,
    "Mus_to_Pek": 0,
    "Mus_to_Mus": 0
}

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

subprocess.run(["samtools", "index", pek_to_mus_bam], check=True)
subprocess.run(["samtools", "index", mus_to_pek_bam], check=True)

print(f"{sample}\t{counts['Pek_to_Pek']}\t{counts['Pek_to_Mus']}\t{counts['Mus_to_Pek']}\t{counts['Mus_to_Mus']}")

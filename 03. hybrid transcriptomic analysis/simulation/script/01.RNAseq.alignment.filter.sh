#!/bin/bash
# RNA-seq alignment pipeline
# Tools: hisat2(v2.1.0), samtools(v1.9)

set -euo pipefail

# default
THREADS=5
OUTDIR=.
SAMPLE="sample"

usage() {
  echo "Usage: $0 -1 read1.fq.gz -2 read2.fq.gz -x genome_index -s sample_name [-t threads]"
  echo ""
  echo "Options:"
  echo "  -1   Read1 fastq.gz"
  echo "  -2   Read2 fastq.gz"
  echo "  -x   Hisat2 genome index prefix"
  echo "  -s   Sample name"
  echo "  -t   Threads (default: 5)"
  exit 1
}

# parse parameters
while getopts "1:2:x:g:a:s:t:h" opt; do
  case $opt in
    1) READ1=$OPTARG ;;
    2) READ2=$OPTARG ;;
    x) GENOME_INDEX=$OPTARG ;;
    s) SAMPLE=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    h) usage ;;
    *) usage ;;
  esac
done

# check required parameters
if [ -z "${READ1:-}" ] || [ -z "${READ2:-}" ] || [ -z "${GENOME_INDEX:-}" ] || [ -z "${SAMPLE:-}" ]; then
  echo "Error: Missing required arguments"
  usage
fi

# set tool path
HISAT2=/gpfs/home/xumingmin/dir.xmm/bin/hisat2-2.1.0/hisat2
SAMTOOLS=/gpfs/CLsoftware/samtools-1.9/bin/samtools

mkdir -p ${SAMPLE}
cd ${SAMPLE}

echo "### Step 1: Hisat2 alignment"
${HISAT2} -p ${THREADS} \
  -x ${GENOME_INDEX} \
  -1 ${READ1} \
  -2 ${READ2} \
  2> ${SAMPLE}.hisat2.log \
  | ${SAMTOOLS} view -bS - > ${SAMPLE}.hisat2.bam

echo "### Step 2: Raw alignment statistics"
${SAMTOOLS} flagstat -@ ${THREADS} ${SAMPLE}.hisat2.bam > ${SAMPLE}.hisat2.raw.flagstat.txt

echo "### Step 3: Filter multiple-mapping and unpaired reads"
${SAMTOOLS} view -q 60 -F 256 -F 12 ${SAMPLE}.hisat2.bam -b > ${SAMPLE}.hisat2.filter.bam

echo "### Step 4: Sort BAM"
${SAMTOOLS} sort -@ ${THREADS} ${SAMPLE}.hisat2.filter.bam -o ${SAMPLE}.hisat2.filter.sorted.bam
rm ${SAMPLE}.hisat2.filter.bam

echo "### Step 5: Filtered alignment statistics"
${SAMTOOLS} flagstat -@ ${THREADS} ${SAMPLE}.hisat2.filter.sorted.bam > ${SAMPLE}.hisat2.filter.flagstat.txt

echo "### Pipeline finished for ${SAMPLE}!"

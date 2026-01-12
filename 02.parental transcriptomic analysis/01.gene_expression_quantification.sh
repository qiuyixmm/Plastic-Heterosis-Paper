#### parental duck species (e.g. Mus.AL.104)
## genome index
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2-build -p 5 Muscovyduck.chromosome_level.genome.fasta Muscovyduck.chromosome_level.genome

## hisat2 alignment
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2 -p 5 \
 -x /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Mus/Muscovyduck.chromosome_level.genome \
 -1 /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/data/Mus.AL/Mus.AL.104_1.fastq.gz \
 -2 /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/data/Mus.AL/Mus.AL.104_2.fastq.gz \
 | /GS01/software/bin/samtools view -bS - > Mus.AL.104.hisat2.bam

## filter multiple mapping reads
## filter unpaired reads
## sort filter bam
/GS01/software/bin/samtools view -q 60 -F 256 -F 12 Mus.AL.104.hisat2.bam -b > Mus.AL.104.hisat2.filter.bam
/GS01/software/bin/samtools sort -@ 5 Mus.AL.104.hisat2.filter.bam -o Mus.AL.104.hisat2.filter.sorted.bam

## reads counts
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/subread/bin/featureCounts -T 5 \
 -a /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Mus/Muscovyduck.genome.modified.gtf \
 -G /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Mus/Muscovyduck.chromosome_level.genome.fasta \
 -J -p -B -C -t exon -g gene_id \
 -o Mus.AL.104.counts.txt Mus.AL.104.hisat2.filter.sorted.bam

## stringtie expression
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/stringtie_2.0.4/stringtie Mus.AL.104.hisat2.filter.sorted.bam \
 -o /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/align/Mus.AL/Mus.AL.104.gtf -p 5 \
 -G /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Mus/Muscovyduck.genome.modified.gtf \
 -A Mus.AL.104.gene_abund.tab -B -e

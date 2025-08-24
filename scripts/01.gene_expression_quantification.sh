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


#### Hybrid duck (e.g. Mul.AL.118)
## hybrid genome index
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2-build -p 5 Pekingduck.Muscovyduck.genome.fa Pekingduck.Muscovyduck.genome

## hisat2 alignment
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2 -p 5 \
 -x /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Pek.Mus/Pekingduck.Muscovyduck.genome \
 -1 /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/data/Mul.AL/Mul.AL.118_1.fastq.gz \
 -2 /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/data/Mul.AL/Mul.AL.118_2.fastq.gz \
 | /GS01/software/bin/samtools view -bS - > Mul.AL.118.hisat2.bam

## filter multiple mapping reads
## filter unpaired reads
## sort filter bam
/GS01/software/bin/samtools view -q 60 -F 256 -F 12 Mul.AL.118.hisat2.bam -b > Mul.AL.118.hisat2.filter.bam
/GS01/software/bin/samtools sort -@ 5 Mul.AL.118.hisat2.filter.bam -o Mul.AL.118.hisat2.filter.sorted.bam

## reads counts
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/subread/bin/featureCounts -T 5 \
 -a /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.gtf \
 -G /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.fa \
 -J -p -B -C -t exon -g gene_id \
 -o Mul.AL.118.counts.txt Mul.AL.118.hisat2.filter.sorted.bam

## stringtie expression
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/stringtie_2.0.4/stringtie Mul.AL.118.hisat2.filter.sorted.bam \
 -o /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/align/Mul.AL/Mul.AL.118.gtf -p 5 \
 -G /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.gtf \
 -A Mul.AL.118.gene_abund.tab -B -e


#### goose
## genome index
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2-build -p 6 goose.assembly.fasta goose.assembly

#### hisat2 alignment
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2 -p 6 \
 -x /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/02.genome/goose.assembly \
 -1 /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/01.data/Gos.AL.69_R1.clean.fastq.gz \
 -2 /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/01.data/Gos.AL.69_R2.clean.fastq.gz \
 | /GS01/software/bin/samtools view -bS | /GS01/software/bin/samtools sort -@ 5 - -o Gos.AL.69.hisat2.sorted.bam

#### reads counts
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/subread/bin/featureCounts -T 5 \
 -a /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/02.genome/goose.modified_annotation.gtf \
 -G /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/02.genome/goose.assembly.fasta \
 -J -p -B -C -t exon -g gene_id \
 -o Goose.AL.OV.counts.txt *.hisat2.sorted.bam

#### stringtie expression
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/stringtie_2.0.4/stringtie Gos.AL.69.hisat2.sorted.bam \
 -o /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/03.align/Gos.AL.69.gtf -p 5 \
 -G /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/02.genome/goose.modified_annotation.gtf \
 -A Gos.AL.69.gene_abund.tab -B -e

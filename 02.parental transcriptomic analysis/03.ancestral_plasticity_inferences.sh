#### goose
## genome index
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/hisat2-2.1.0/hisat2-build -p 6 goose.assembly.fasta goose.assembly

#### hisat2 alignment (eg. Gos.AL.69)
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


#### Multiple branch estimates
## For ab libitum feeding condition
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/python3.9/envs/CAGEE/bin/cagee --cores 5 \
 --tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/tree.txt \
 --sigma_tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/tree.sigma.alter.txt \
 --infile /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/Pek.Mus.Gos.AL.expression.TPM.txt \
 --output_prefix CAGEE.res.alter

## For overfeeding condition
/GS01/project/pengms_group/pengms20t1/dir.xumm/bin/python3.9/envs/CAGEE/bin/cagee --cores 5 \
 --tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/02.OV/tree.txt \
 --sigma_tree /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/02.OV/tree.sigma.alter.txt \
 --infile /GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/goose/05.CAGEE/02.OV/Pek.Mus.Gos.OV.expression.TPM.txt \
 --output_prefix CAGEE.res.alter

cd /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/01.AL/05.read_assign

## for parental species
/gpfs/home/xumingmin/dir.xmm/bin/subread-1.6.2/bin/featureCounts -T 5 -R CORE \
 -a /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.gtf \
 -G /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.fa \
 -J -p -B -C -t exon -g gene_id -o parental.AL.counts.txt \
 /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/align/Pek.AL/*.hisat2.filter.sorted.bam \
 /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/align/Mus.AL/*.hisat2.filter.sorted.bam 

## For pseudo-hybrids
/gpfs/home/xumingmin/dir.xmm/bin/subread-1.6.2/bin/featureCounts -T 5 -R CORE \
 -a /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.gtf \
 -G /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.fa \
 -J -p -B -C -t exon -g gene_id -o Hybrids.AL.counts.txt \
 /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/01.AL/02.align/*/*.hisat2.filter.sorted.bam


cd /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/02.OV/05.read_assign

## for parental species
/gpfs/home/xumingmin/dir.xmm/bin/subread-1.6.2/bin/featureCounts -T 5 -R CORE \
 -a /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.gtf \
 -G /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.fa \
 -J -p -B -C -t exon -g gene_id -o parental.AL.counts.txt \
 /gpfs/home/xumingmin/dir.xmm/mule_duck/overfed/align/Pek.OV/*.hisat2.filter.sorted.bam \
 /gpfs/home/xumingmin/dir.xmm/mule_duck/overfed/align/Mus.OV/*.hisat2.filter.sorted.bam

## For pseudo-hybrids
/gpfs/home/xumingmin/dir.xmm/bin/subread-1.6.2/bin/featureCounts -T 5 -R CORE \
 -a /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.gtf \
 -G /gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome.fa \
 -J -p -B -C -t exon -g gene_id -o Hybrids.AL.counts.txt \
 /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/02.OV/02.align/*/*.hisat2.filter.sorted.bam


## summary
cd /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/01.AL/05.read_assign
bash AL_make_summary.sh

cd /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/02.OV/05.read_assign
bash OV_make_summary.sh

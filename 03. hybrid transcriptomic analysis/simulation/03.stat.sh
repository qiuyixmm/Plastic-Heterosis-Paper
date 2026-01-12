## for ad libitum feeding conditions
for bam in /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/01.AL/02.align/*/*.sorted.bam; do
    /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/tools/03.count_alignment_origin.py $bam >> AL.alignment_summary.tsv
done


## for overfeeding conditions
for bam in /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/02.OV/02.align/*/*.sorted.bam; do
    /gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/tools/03.count_alignment_origin.py $bam >> OV.alignment_summary.tsv
done

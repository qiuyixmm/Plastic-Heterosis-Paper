# path to genome index
genome="/gpfs/home/xumingmin/dir.xmm/mule_duck/control.fed/genome/Pek.Mus/Pekingduck.Muscovyduck.genome"

# path to script
script="/gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/tools/01.RNAseq.alignment.filter.sh"


## for ad libitum feeding condition
# sample list
samples=(
  Pek.AL.101.Mus.AL.96
  Pek.AL.105.Mus.AL.88
  Pek.AL.109.Mus.AL.112
  Pek.AL.113.Mus.AL.116
  Pek.AL.117.Mus.AL.120
  Pek.AL.81.Mus.AL.108
  Pek.AL.85.Mus.AL.84
  Pek.AL.93.Mus.AL.92
  Pek.AL.97.Mus.AL.104
)


# data directory
data_dir="/gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/01.AL/01.data"

# run
for sample in "${samples[@]}"; do
    echo "Running sample ${sample}..."
    ${script} -t 5 \
      -x ${genome} \
      -1 ${data_dir}/${sample}_1.fastq.gz \
      -2 ${data_dir}/${sample}_2.fastq.gz \
      -s ${sample}
    echo "Sample ${sample} finished."
    echo "--------------------------------------------"
done


## for overfeeding condition
# sample list
samples=(
  Pek.OV.13.Mus.OV.16
  Pek.OV.17.Mus.OV.20
  Pek.OV.1.Mus.OV.32
  Pek.OV.21.Mus.OV.4
  Pek.OV.25.Mus.OV.36
  Pek.OV.29.Mus.OV.8
  Pek.OV.33.Mus.OV.40
  Pek.OV.37.Mus.OV.28
  Pek.OV.5.Mus.OV.24
  Pek.OV.9.Mus.OV.12
)

# data directory
data_dir="/gpfs/home/xumingmin/dir.xmm/mule_duck/01.major_revison/04.simulation/02.OV/01.data"

# run
for sample in "${samples[@]}"; do
    echo "Running sample ${sample}..."
    ${script} -t 5 \
      -x ${genome} \
      -1 ${data_dir}/${sample}_1.fastq.gz \
      -2 ${data_dir}/${sample}_2.fastq.gz \
      -s ${sample}
    echo "Sample ${sample} finished."
    echo "--------------------------------------------"
done

#$ -t 1-6
#$ -q lupiengroup
#$ -V
#$ -cwd
#$ -N ReadsInConsensus

BAMS=( \
    "../../Data/Processed/2018-09-17/output_Ctrl/align/rep1/Rashimi_ATAC1_S25_L007_R1_001.nodup.bam" \
    "../../Data/Processed/2018-09-17/output_Ctrl/align/rep2/Rashim_ATAC2_S26_L007_R1_001.nodup.bam" \
    "../../Data/Processed/2018-09-17/output_Ctrl/align/rep3/Rashim_ATAC3_S27_L007_R1_001.nodup.bam" \
    "../../Data/Processed/2018-09-17/output_MB6/align/rep1/Rashim_ATAC4_S28_L007_R1_001.nodup.bam" \
    "../../Data/Processed/2018-09-17/output_MB6/align/rep2/5_S21_L006_R1_001.nodup.bam" \
    "../../Data/Processed/2018-09-17/output_MB6/align/rep3/6_S22_L006_R1_001.nodup.bam"
)
CONDITIONS=("Ctrl" "Ctrl" "Ctrl" "MB6" "MB6" "MB6")
REPS=(1 2 3 1 2 3)

i=$(($SGE_TASK_ID - 1))

# ensure bedtools v 2.23.0
bedtools --version

con=${CONDITIONS[$i]}
rep=${REPS[$i]}
bam=${BAMS[$i]}
echo "$con Rep$rep"
consensus="../2018-11-28_global-accessbility/Consensus/consensus.Ctrl-MB6.bed"

# calculate read depth in each peak of the consensus peak list
echo "Calculating depths"
echo "bedtools coverage -abam $bam -b $consensus -counts > Counts/${con}_Rep${rep}.bed"
bedtools coverage -abam $bam -b $consensus -counts > Counts/${con}_Rep${rep}.bed

echo "Sorting"
echo "sort -k1,1 -k2,2n Counts/${con}_Rep${rep}.bed > Counts/${con}_Rep${rep}.sorted.bed"
sort -k1,1 -k2,2n Counts/${con}_Rep${rep}.bed > Counts/${con}_Rep${rep}.sorted.bed
echo "Done"

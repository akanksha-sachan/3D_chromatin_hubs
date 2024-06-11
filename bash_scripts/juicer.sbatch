#!/bin/bash
#SBATCH -p ma2-gpu
#SBATCH --job-name juicer_tools
#SBATCH --output juicer-log-%J.txt
#SBATCH -w compute-4-13
#SBATCH --time=1-00:00:00

#load cuda
module load cuda-12.2

#activate env which has java
source /home/asachan2/miniconda3/etc/profile.d/conda.sh
conda activate hubs

#get pc1 for specific res and chrom, -p flag necessary for resolutions>500kb
res=100000
res0=100kb
chrom=1

java -jar /work/magroup/asachan2/genomic_hubs/juicer_tools/juicer_tools_1.19.02.jar eigenvector -p VC /work/magroup/asachan2/genomic_hubs/data/GM12878/4DNData/4DNFI9YAVTI1.hic ${chrom} BP ${res} /work/magroup/asachan2/genomic_hubs/tmp/GM12878/pc_one_${res0}.txt


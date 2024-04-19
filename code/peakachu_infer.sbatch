#!/bin/bash
#SBATCH -p ma2-gpu
#SBATCH --job-name peakachu
#SBATCH --output job-log-%J.txt
#SBATCH -w compute-4-26
#SBATCH --time=1-00:00:00

#Load CUDA
module load cuda-12.2

# --> uncomment to use conda environment with jupyter installed <---
source /home/asachan2/miniconda3/etc/profile.d/conda.sh
conda activate hubs

python score_chromosome.py --resolution 10000 --chrom chr2 --model /work/magroup/asachan2/genomic_hubs/peakachu/models/GM12878/chr2.pkl --path /work/magroup/asachan2/genomic_hubs/data/GM12878/4DNFITRVKRPA_GM12878.mcool::/resolutions/10000 --output /work/magroup/asachan2/genomic_hubs/peakachu/predicted_loops/GM12878-H3K27ac-chr2-scores.bedpe

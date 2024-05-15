#!/bin/bash
#SBATCH -p ma2-gpu
#SBATCH --job-name jupyter_notebook
#SBATCH --output jupyter-log-%J.txt
#SBATCH -w compute-4-26
#SBATCH --time=1-00:00:00

#Load CUDA
module load cuda-12.2

# --> uncomment to use conda environment with jupyter installed <---
source /home/asachan2/miniconda3/etc/profile.d/conda.sh
conda activate hubs

###### RUNNING JUPYTERLAB ########
python -m ipykernel install --user --name hubs --display-name "Python (hubs)"
## get tunneling info
XDG_RUNTIME_DIR=""
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)

## print tunneling instructions to sbatch.stdout
echo -e "
    Copy/Paste this in your local terminal to ssh tunnel with remote
    -----------------------------------------------------------------
    ssh -N -L $ipnport:$ipnip:$ipnport asachan2@lanec1.compbio.cs.cmu.edu
    -----------------------------------------------------------------
    Then open a browser on your local machine to the following address
    ------------------------------------------------------------------
    localhost:$ipnport  (prefix w/ https:// if using password)
    ------------------------------------------------------------------
    "


# start an ipcluster instance and launch jupyter server
jupyter-lab --NotebookApp.iopub_data_rate_limit=100000000000000 --no-browser --port=$ipnport --ip=$ipnip --NotebookApp.password='' --NotebookApp.token='' --notebook-dir=/work/magroup/asachan2/genomic_hubs

#!/bin/bash
#SBATCH --job-name code-server
#SBATCH --time 0-12:00:00
#SBATCH --output=vscode-%J.log 

#SBATCH --mem=30Gb
#SBATCH -p ma2-gpu
#SBATCH -w compute-4-26
#SBATCH --ntasks-per-node=1

 
PASSWORD=29c92b95d5a121a897379cd1 
PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)

echo "********************************************************************" 
echo "Starting code-server in Slurm"
echo "Environment information:" 
echo "Date:" $(date)
echo "Allocated node:" $(hostname)
echo "Node IP:" $(hostname -i)
echo "Path:" $(pwd)
echo "Password to access VSCode:" $PASSWORD
echo "********************************************************************" 

echo -e "
    Copy/Paste this in your local terminal to ssh tunnel with remote
    -----------------------------------------------------------------
    ssh -N -L $ipnport:$ipnip:$ipnport $(whoami)@lanec1.compbio.cs.cmu.edu
    -----------------------------------------------------------------
    Then open a browser on your local machine to the following address
    ------------------------------------------------------------------
    localhost:$ipnport  (prefix w/ https:// if using password)
    ------------------------------------------------------------------
    "
PASSWORD=$PASSWORD ~/.local/lib/code-server-4.7.0/bin/code-server --bind-addr 0.0.0.0:$ipnport --auth password --disable-telemetry

#code to assign bins to anchors etc, other things to be done to all feature calls etc

#need a function to make all the temp directories per cell line

#function for creating bins on the reference genome using bedtools python wrapper

import subprocess
import os
import pandas as pd
import numpy as np

def bedtools_makewindows(chromsizes_file, resolution, tmp_dir):
    """runs this command: $ bedtools makewindows -g hg19.txt -w 1000000
        can add additional if statements to be able to run -b input.bed and -n 10 command
    """

    #create sub-directory for resoltuion/window sizes 
    resolution_dir = os.path.join(tmp_dir, str(resolution))
    if not os.path.exists(resolution_dir):
        os.makedirs(resolution_dir)

    #extract name of output file from chromsizes_file (example: bins_hg19.bed)
    output_file = f'bins_{os.path.basename(chromsizes_file).split(".")[0]}.bed'
    output_path = os.path.join(resolution_dir, output_file)

    # Construct the bedtools command
    cmd = f'bedtools makewindows -g "{chromsizes_file}" -w {resolution}'

    # Run the command and redirect output to a file
    with open(output_path, 'w') as file:
        process = subprocess.run(cmd, shell=True, stdout=file, stderr=subprocess.PIPE, text=True)

    # Handle command output
    if process.returncode == 0:
        print(f"Successfully created bins file: {output_path}")
    else:
        print(f"Error creating bins file: {output_file}. Error: {process.stderr}")

def convertGenCoordsToBinPos(genCoords, resolution):
    """
    Convert genomic coordinates (base pairs) to bin positions (indices used in matrix representation) of a given resolution
    Assumption: genomic coordinates are from the same chromosome
    """
    regionIndices = [0]*len(genCoords)
    for i in range(len(genCoords)):
        regionIndices[i] = genCoords[i]//resolution
    return regionIndices


#test function
chromsizes_file = '/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/hg38.chrom.sizes'
tmp_dir = '/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp'
resolution = 10000
bedtools_makewindows(chromsizes_file, resolution, tmp_dir)
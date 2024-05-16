######### GENE EXPRESSION TO CHROMATIN CONTACT NODE MAPPING #########

# input: GTF file to locate gene TSS on the human genome, RNA-seq data, cooler of specified resolution
# output: {gene_name: bin_ID in the cooler at a specific resolution}
# creating intermediate files is okay here as this processing is only done once, at the start of the pipeline; make objects and use class methods to manipulate other data repeatedly instead of using tmp files

import pandas as pd
import cooler 
import numpy as np
import os
import ipdb
import requests, sys
import json
from tqdm import tqdm

def fetch_gene_ensemblAPI(ensembl_ids, batch_size=1000, ref_genome='hg38'):
    #add 'chr' in front of chrom number for convention standard whenever using ensemble 
    """
    Fetch details for multiple Ensembl IDs using api 
    Parameters:
        ensembl_ids (list): A list of Ensembl gene IDs (e.g., ENSG00000001497)
        batch_size (int): The number of IDs to query in each batch
    Returns:
        DataFrame: A DataFrame containing 'id', 'biotype', 'chrom', 'Txstart', 'Txend', 'gene_name' for each gene
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    gene_details = [] 

    #process batches of Ensembl IDs
    for i in tqdm(range(0, len(ensembl_ids), batch_size)):
        batch_ids = ensembl_ids[i:i+batch_size]  #current batch
        data = json.dumps({"ids": batch_ids})

        response = requests.post(server + ext, headers=headers, data=data)
        if response.ok:
            results = response.json()
            #extract details for each gene in the batch
            batch_details = [{
                'id': result['id'],
                'biotype': result.get('biotype'),
                'chrom': result.get('seq_region_name'),
                'Txstart': result.get('start'),
                'Txend': result.get('end'),
                'gene_name': result.get('display_name')
            } for result in results.values() if 'id' in result]
            gene_details.extend(batch_details) 
        else:
            print(f"API Request Failed: {response.status_code} {response.reason} for batch starting at index {i}")

    df = pd.DataFrame(gene_details)  
    return df

def fetch_gene_gtf(gtf_file, ensembl_ids, ref_genome='hg38'):
    #the headers are named from GTF convention: universal function for all total RNA-seq data for hg38
    """
    Fetch gene details from a GTF file for a list of Ensembl IDs
    Parameters:
        gtf_file (str): Path to the GTF file
        ensembl_ids (list): A list of Ensembl gene IDs (e.g., ENSG00000001497); 
        we don't input the version of the gene annotation (post decimal part of the ensembl id) as the gtf file can vary for different datasets
        assumption is that the gene details needed for this project (tss start and end) doen't vary much between gene annotation versions (only exons etc are updated)
    Returns:
        pandas df: A DataFrame containing 'gene_id', 'gene_type', 'chrom', 'Txstart', 'Txend', 'gene_name' for each gene
    """
    #load gtf file
    gtf = pd.read_csv(gtf_file, sep='\t', header=None, comment='#') #skips the commented rows
    gtf_genes = gtf[gtf[2] == 'gene'].copy() #only select rows with gene annotation
    gtf_genes['gene_id'] = gtf_genes[8].str.extract(r'gene_id "([^"]+)"') 
    gtf_genes['gene_name'] = gtf_genes[8].str.extract(r'gene_name "([^"]+)"') 
    gtf_genes['gene_type'] = gtf_genes[8].str.extract(r'gene_type "([^"]+)"') 
    #assign columns as chrom, Txstart, Txend headers
    gtf_genes['chrom'] = gtf_genes[0]
    gtf_genes['Txstart'] = gtf_genes[3]
    gtf_genes['Txend'] = gtf_genes[4]
    df = gtf_genes[['gene_id', 'gene_type', 'gene_name', 'chrom', 'Txstart', 'Txend']].copy()

    #select only those genes that are in the ensembl_ids list
    #strip the gene_id at . to query the ensembl_ids in the df // this enables us to query ids from any version of the gtf
    df['gene_id'] = df['gene_id'].str.split('.').str[0]
    fetched_genes_df = df[df['gene_id'].isin(ensembl_ids)]
    return fetched_genes_df

def processRNAseqData(infile, outfile, gtf_file, threshold=0.1, normalisation='log1p_t', ref_genome='hg38'):
    #the headers are custom to the current rna-seq data: specific function, switch headers according to input file
    """
    Read in raw total RNA-seq file for gene quantification, select a normalisation method and return the normalised data as gene_name and normalised and thresholded gene expression
    Parameters:
        infile (str): Path to the raw RNA-seq file
    Returns:
        file csv; column headers: gene names and expressed gene expression values
    """
    #load total RNA-seq data, with ensemble ids as a column, strip the ids at . and make a list of id strings to pass to the API
    rna_data = pd.read_csv(infile, sep='\t')
    #only select columns with gene_id starting from ENSG and expression values (pme_TPM)
    rna_data = rna_data[['gene_id', 'pme_TPM']]
    #only select ensemble ids and strip them before decimal point for futher processing
    rna_data = rna_data[rna_data['gene_id'].str.startswith('ENSG')]
    rna_data['gene_id'] = rna_data['gene_id'].str.split('.').str[0]
    #drop duplicates
    rna_data = rna_data.drop_duplicates(subset='gene_id')
    #ensure pme_TPM is float
    rna_data['pme_TPM'] = rna_data['pme_TPM'].astype(float)
    ensembl_ids = rna_data['gene_id'].tolist()

    #query the ensembl ids of the rna expression data in the gtf file for additional details
    gene_df = fetch_gene_gtf(gtf_file, ensembl_ids)
    #only select those genes that have gene_type as protein_coding
    gene_df = gene_df[gene_df['gene_type'] == 'protein_coding']
    #drop duplicates
    gene_df = gene_df.drop_duplicates(subset='gene_id')

    #merge gene_df with rna expression data on ensembl_id to keep info about protein coding genes
    expressed_genes = pd.merge(rna_data, gene_df, left_on='gene_id', right_on='gene_id', how='inner')

    #normalise the expression values, can include additional normalisations
    if normalisation == 'log1p_t':
        expressed_genes['normalized_expression'] = np.log1p(expressed_genes['pme_TPM'])
    else:
        expressed_genes['normalized_expression'] = expressed_genes['pme_TPM']

    #TODO: threshold the expression values based on differential expression?

    #sort rows by chr string (1 is top)
    expressed_genes = expressed_genes.sort_values(by='chrom')
    expressed_genes.to_csv(outfile, sep='\t', index=False) #saves df without the index
    return None #performs only I/O

def assignBins2Genes(infile, bins_files, outfile, ref_genome='hg38'):
    """
    Read in the expressed gene list and assign chromosome locations to the gene TSSs, and gene body
    Parameters:
        infile (str): Path to the file containing gene list of expressed genes
        bin_files (dict): list of bin bed files (values) at multiple resolutions (keys) to get chr:start-end for each bin
        ref_genome (str): Reference genome to use (e.g., 'hg19', 'hg38', etc.)
    Returns:
        file csv; column headers: gene_name, TSSBin_10kb: chrN:start, TSSBin_{other}, BodyBins_10kb: tuple of chrN:starts, BodyBins_{other}: tuple 
    """
    #load the gene_name, and start end locations
    genes = pd.read_csv(infile, sep='\t')
    #copy gene_name chrom txStart txEnd normalized_expression to output df
    gene_chrom_bin_num = genes[['gene_name', 'chrom', 'Txstart', 'Txend', 'normalized_expression']].copy()

    for res0, bins_file in bins_files.items():
        #load the bins file
        bins = pd.read_csv(bins_file, sep='\t', header=None, names=['chrom', 'start', 'end'])

        #dict for intervalIndices for each chromosome
        chrom_intervals = {}
        for chrom, group in bins.groupby('chrom'):
            #create intervals for all bins per chromosome; starts belong to the interval, ends are excluded
            chrom_intervals[chrom] = pd.IntervalIndex.from_arrays(group['start'], group['end'], closed='left')
        
        #assign a bin to the TSS per gene; input a row from the gene_list
        def getTSSBin(row):
            if row['chrom'] in chrom_intervals:
                intervals = chrom_intervals[row['chrom']]
                #check if there's a direct interval containing the TSS
                try:
                    idx = intervals.get_loc(row['Txstart'])
                    interval = intervals[idx]
                    return f"{row['chrom']}:{interval.left}"
                except KeyError:
                    #if no direct interval, find the closest preceding interval (manual 'ffill')
                    preceding_intervals = intervals[intervals.left <= row['Txstart']]
                    if not preceding_intervals.empty:
                        interval = preceding_intervals[-1]  # Last interval before or at TSS
                        return f"{row['chrom']}:{interval.left}"
            return None
        gene_chrom_bin_num[f'TSSBin_{res0}'] = gene_chrom_bin_num.apply(getTSSBin, axis=1)

        #assign bins to the gene body
        def getBodyBins(row):
            if row['chrom'] in chrom_intervals:
                intervals = chrom_intervals[row['chrom']]
                overlapping_intervals = intervals.overlaps(pd.Interval(row['Txstart'], row['Txend'], closed='right'))
                return tuple(f"{row['chrom']}:{intv.left}" for intv in intervals[overlapping_intervals])
            return ()
        gene_chrom_bin_num[f'BodyBins_{res0}'] = gene_chrom_bin_num.apply(getBodyBins, axis=1)

    #save the bin location for each gene to a file
    gene_chrom_bin_num.to_csv(outfile, sep='\t', index=False)
    return None #performs only I/O
    
def assignTSS2HiCBinID(infile, res0, bins_3D, ref_genome='hg38'):
    #make sure the ref genomes for both rna and hic are the same
    #TODO: debug for new gene list type and dataset
    """
    Reads a list of genes with bin starts defined at different resolutions and assigns the gene to bin_ID of the current resolution at the whole genome level
    Parameters:
        infile (str): Path to the file containing the gene list.
        resolution (str): The resolution key to use (e.g., '1kb', '5kb', etc.)
        bins (cooler dataFrame): bin_id in hic and its genomic location
    Returns:
        dict: A dictionary mapping each gene to its bin identifier (int) to index the cooler matrix at a specific resolution
    """
    # Read the gene list file, it has headers like Bin_10k, Bin_5k, etc
    gene_list = pd.read_csv(infile, sep='\t')
    TSSbin_res = f"TSSBin_{res0}" #header of infile containing bin info of the TSS

    # Make gene_bins df
    gene_bins = gene_list[['gene_name', TSSbin_res]].copy()
    # Split 'bin_res' into 'chrom' and 'start', ensuring correct assignment 
    gene_bins[['chrom', 'start']] = gene_bins[TSSbin_res].str.split(':', expand=True)
    gene_bins['start'] = gene_bins['start'].astype(str).str.replace(',', '').astype(int)  # Convert 'start' to integer

    # Name index column of bins df as bin_id
    bins_3D = bins_3D.reset_index()
    bins_3D = bins_3D.rename(columns={'index': 'bin_id'})
    # Select bin_id, chrom, start columns and ensure start is integer
    bins_3D['start'] = bins_3D['start'].astype(int)
    
    # Merge the DataFrames on 'chrom' and 'start' using a left join
    merged_df = pd.merge(gene_bins, bins_3D, on=['chrom', 'start'], how='left')

    # Assign NaN to genes that are out of bounds for cooler bins (chromsizes aren't fully captured in cooler bins)
    non_matched_genes = merged_df[merged_df['bin_id'].isnull()]
    
    # Convert to dictionary
    gene_to_bin_dict = dict(zip(merged_df.loc[~merged_df['gene_name'].isin(non_matched_genes['gene_name']), 'gene_name'],
                            merged_df.loc[~merged_df['gene_name'].isin(non_matched_genes['gene_name']), 'bin_id']))
    return gene_to_bin_dict , non_matched_genes 

# multi-res assignment of genes to bins (gene to binID file already has multi-res bin_ids that can be assigned if the cooler bins were multi-res)


    
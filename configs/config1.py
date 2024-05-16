#config file for genomic hubs workflow

class GenomicParams:
    def __init__(self):
        self.cell_types = ["GM12878", "K562"] 
        self.ref_genome = "hg38"
        self.resolutions = [10000, 1000000]
        self.res_strs = ["10kb", "1Mb"]
        self.chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

class Paths:
    def __init__(self, cell_type, resolution, ref_genome):
        self.gtf_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/gencode.v45.basic.annotation.{ref_genome}.gtf"
        self.rnaseq_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{cell_type}/ENCODE_hg38/RNA_ENCSR820PHH/ENCFF345SHY_gene_quant.tsv"
        self.expressed_genes_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/expressed_genes_{ref_genome}_{cell_type}.txt"
        self.gene_chrom_bin_num_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/gene_chrom_bin_num_{ref_genome}_{cell_type}.txt"
        self.cool_file = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{cell_type}/4DN Data Portal_hg19/4DNFITRVKRPA_GM12878.mcool"
        self.hic_file = "https://www.encodeproject.org/files/ENCFF318GOM/@@download/ENCFF318GOM.hic"
        self.bins_file = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{resolution}/bins_{ref_genome}.bed"
        self.temp_dir = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/code/temp/{cell_type}"

class Config:
    def __init__(self, cell_type, resolution):
        self.genomic_params = GenomicParams()
        self.paths = Paths(cell_type, resolution, self.genomic_params.ref_genome)

if __name__ == "__main__":
    config = Config("GM12878", 10000)
    print(config.paths.gtf_infile)
    print(config.paths.gene_chrom_bin_num_outfile)
    print(config.paths.bins_file)

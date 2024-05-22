# use scratch folder for loading data on the cluster

class Config:
    """Main config class to store parameters and paths for the workflow"""

    def __init__(self, cell_type, resolution):
        self.genomic_params = self.genomic_params()
        self.paths = self.paths(cell_type, resolution, self.genomic_params.ref_genome)

    class genomic_params:
        """Nested class to store genomic parameters for the workflow"""

        def __init__(self):
            self.cell_types = ["GM12878"]
            self.ref_genome = "hg38"  # used to read paths dynamically in config file
            self.resolutions = [1000000, 10000]  # list of integers to query
            self.res_strs = [
                "1Mb",
                "10kb",
            ]  # list of strings based on standard naming convention of resolutions
            self.hic_norm = "VC"  # vanilla coverage normalization
            self.chromosomes = [
                "chr1",
                "chr2",
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chrX",
                "chrY",
            ]  # list of strings based on standard naming convention of chromosomes

    class paths:
        """Nested class to store paths for the workflow"""

        def __init__(self, cell_type, resolution, ref_genome):
            self.cell_type = cell_type
            self.resolution = resolution
            self.ref_genome = ref_genome
            self.initialize_rnaseq_paths()
            self.initialize_bed_paths()
            self.initialize_hic_paths()
            self.initialize_temp_paths()

        def initialize_rnaseq_paths(self):
            """Initialize paths for RNA-seq data"""
            self.gtf_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/gencode.v45.basic.annotation.{self.ref_genome}.gtf"
            self.rnaseq_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/ENCODE_4/RNA_ENCSR820PHH/ENCFF345SHY_gene_quant.tsv"
            self.expressed_genes_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/expressed_genes_{self.ref_genome}_{self.cell_type}.txt"
            self.gene_chrom_bin_num_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/gene_chrom_bin_num_{self.ref_genome}_{self.cell_type}.txt"

        def initialize_bed_paths(self):
            """Initialize paths for indexing genomic sequences"""
            self.bins_file = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.resolution}/bins_{self.ref_genome}.bed"
            self.chrom_sizes_file = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.ref_genome}.chrom.sizes"
        
        def initialize_hic_paths(self):
            """Initialize paths for Hi-C data"""
            self.cool_file = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DN_dataportal_{self.ref_genome}/4DNFITRVKRPA.mcool"
            self.hic_file = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DN_dataportal_{self.ref_genome}/4DNFI9YAVTI1.hic"

        def initialize_temp_paths(self):
            """Initialize paths for temporary files"""
            self.temp_dir = (
                f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.cell_type}"
            )

"""config file to specify params for runs"""


class Config:
    """Main config class to store parameters and paths for the workflow"""

    def __init__(self):
        self.genomic_params = self.genomic_params()
        self.paths = self.paths(self.genomic_params.cell_types[0])
    
    class genomic_params:
        """Nested class to store genomic parameters for the workflow"""

        def __init__(self):
            self.cell_types = ["GM12878"]
            self.ref_genome = "hg38"  # used to read paths dynamically in config file
            self.resolutions = [1000000, 1000]  # list of integers to query
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
        """Nested class to store input file paths"""

        def __init__(self, cell_type):
            self.cell_type = cell_type
            self.initialize_rnaseq_paths()
            self.initialize_bed_paths()
            self.initialize_hic_paths()
            self.initialize_tool_paths()
            self.initialize_temp_paths()
            self.initialize_loops_paths()
            self.initialize_hub_paths()
            self.initialize_bigwig_paths()

        def initialize_rnaseq_paths(self):
            """Initialize paths for RNA-seq data"""
            self.gtf_infile = "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/gencode.v45.basic.annotation.hg38.gtf"
            self.rnaseq_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/ENCODE_4/RNA_ENCSR820PHH/ENCFF345SHY_gene_quant.tsv"
            self.expressed_genes_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/expressed_genes_hg38_{self.cell_type}.txt"
            self.gene_chrom_bin_num_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/gene_chrom_bin_num_hg38_{self.cell_type}.txt"

        def initialize_bed_paths(self):
            """Initialize paths for indexing genomic sequences"""
            self.ref_genome_bins_infile = "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/10kb/bins_hg38.bed"
            self.chrom_sizes_infile = "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/hg38.chrom.sizes"
            
        def initialize_hic_paths(self):
            """Initialize paths for Hi-C data"""
            self.cool_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DNData/4DNFITRVKRPA.mcool"
            self.hic_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DNData/4DNFI9YAVTI1.hic"
        
        def initialize_loops_paths(self):
            """Initialize paths for Hi-C data inputs and edgelist outputs in bedpe format"""
            self.hiccups_merged_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.cell_type}/loop_hiccups/merged_loops.bedpe"
            self.loops_txt_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DNData/GSE63525_replicate_hg19_HiCCUPS_looplist.txt"
            self.loops_bedpe_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.cell_type}/GSE63525_replicate_hg19_HiCCUPS_looplist.bedpe"
            self.geo_loops_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DNData/GSE63525_replicate_hg38liftOver_HiCCUPS_10kb.bedpe"
        
        def initialize_bigwig_paths(self):
            """Initialize paths for bigwigs"""
            self.compartments_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DNData/4DNFIDPK7WFE_compartments.bw"

        def initialize_tool_paths(self):
            """Initialize paths for other tools used in pipeline"""
            self.juicer_tools = "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/juicer_tools/juicer_tools_1.19.02.jar"

        def initialize_temp_paths(self):
            """Initialize paths for temporary files"""
            self.temp_dir = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.cell_type}"

        def initialize_hub_paths(self):
            """ nodesets of different resolutions stored in hdf5 format """
            self.edgelist_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.cell_type}/edgelist.h5"
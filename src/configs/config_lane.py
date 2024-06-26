"""config file to specify params for runs"""
#set cell type and resolution to load files from 

class Config:
    """Main config class to store parameters and paths for the workflow"""

    def __init__(self):
        self.genomic_params = self.genomic_params()
        #set currents for runs
        self.current_cell_type = self.genomic_params.cell_types[0]
        self.current_res_str = self.genomic_params.res_strs[1]
        self.current_res = self.genomic_params.resolutions[1]
        self.paths = self.paths(self.current_cell_type, self.current_res_str)
    
    class genomic_params:
        """Nested class to store genomic parameters for the workflow"""

        def __init__(self):
            self.cell_types = ["GM12878", "K562"]
            self.ref_genome = "hg38"  # used to read paths dynamically in config file
            self.resolutions = [1000000, 100000, 10000]  # list of integers to query
            self.res_strs = [
                "1Mb",
                "100kb",
                "10kb"
            ]  # list of strings based on standard naming convention of resolutions
            self.affinity_key = f"OE_{self.res_strs[1]}_affinity" #for saving the affinity plots to this folder
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
            self.oe_matrix_viz()
        
        def oe_matrix_viz(self):
            self.threshold = 0.9 # threshold for observed over expected matrix query of edges; tweak based on single chr viz
            self.start = 0 # start index for OE numpy slice for visualization
            self.end = 30000000 # end index for OE numpy slice for visualization

    class paths:
        """Nested class to store input file paths"""

        def __init__(self, cell_type, res_str):
            self.cell_type = cell_type
            self.res_str = res_str
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
            self.gtf_infile = "/work/magroup/asachan2/genomic_hubs/data/gencode.v45.basic.annotation.hg38.gtf"
            self.rnaseq_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/ENCODE/RNA_ENCSR820PHH/ENCFF345SHY_gene_quant.tsv"
            self.expressed_genes_outfile = f"/work/magroup/asachan2/genomic_hubs/tmp/expressed_genes_hg38_{self.cell_type}.txt"
            self.gene_chrom_bin_num_outfile = f"/work/magroup/asachan2/genomic_hubs/tmp/gene_chrom_bin_num_hg38_{self.cell_type}.txt"

        def initialize_bed_paths(self):
            """Initialize paths for indexing genomic sequences"""
            self.ref_genome_bins = f"/work/magroup/asachan2/genomic_hubs/tmp/bins/bins_{self.res_str}_hg38.bed"
            self.chrom_sizes_infile = "/work/magroup/asachan2/genomic_hubs/data/hg38.chrom.sizes"
            self.sub_compartments_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/4DNData/GSE63525_subcompartments_hg38.bed"
            
        def initialize_hic_paths(self):
            """Initialize paths for Hi-C data"""
            self.cool_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/4DNData/4DNFITRVKRPA.mcool"
            self.hic_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/4DNData/4DNFI9YAVTI1.hic"
        
        def initialize_loops_paths(self):
            """Initialize paths for Hi-C data inputs and edgelist outputs in bedpe format"""
            self.hiccups_merged_infile = "/work/magroup/asachan2/genomic_hubs/juicer_tools/results/whole_genome_gm12878/merged_loops.bedpe"
            self.loops_txt_infile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/{self.cell_type}/4DNData/GSE63525_replicate_hg19_HiCCUPS_looplist.txt"
            self.loops_bedpe_outfile = f"/Users/Akanksha/MaGroup/Genomic Hubs/workflow/tmp/{self.cell_type}/GSE63525_replicate_hg19_HiCCUPS_looplist.bedpe"
            self.replicate_lifted_loops_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/4DNData/GSE63525_replicate_hg38liftOver_HiCCUPS_10kb.bedpe"
            self.all_lifted_loops_infile = "/work/magroup/4DN/yang/Rao_et_al/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_hg38.txt"
        
        def initialize_bigwig_paths(self):
            """Initialize paths for bigwigs"""
            self.compartments_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/4DNData/4DNFIDPK7WFE_compartments.bw"
            self.insulation_infile = f"/work/magroup/asachan2/genomic_hubs/data/{self.cell_type}/4DNData/4DNFI62JTGEX_insulation.bw" #cooltools created file

        def initialize_tool_paths(self):
            """Initialize paths for other tools used in pipeline"""
            self.juicer_tools = "/work/magroup/asachan2/genomic_hubs/juicer_tools/juicer_tools_1.19.02.jar"

        def initialize_temp_paths(self):
            """Initialize paths for temporary files"""
            self.temp_dir = f"/work/magroup/asachan2/genomic_hubs/tmp/{self.cell_type}"
            self.gexf_dir = f"/work/magroup/asachan2/genomic_hubs/tmp/{self.cell_type}/graphs/{self.res_str}"

        def initialize_hub_paths(self):
            """ edges of different resolutions stored in hdf5 format """
            self.edgelist_outfile = f"/work/magroup/asachan2/genomic_hubs/tmp/{self.cell_type}/edgelist.h5"
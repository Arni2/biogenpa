import pysam

# File paths
fasta_file = "...plus_decoy_hla.fa"
file_path_i = "...t_plus_decoy_hla.fa.fai"

def get_sequence_around(chromosome, position):
    try:
        # Open the FASTA file
        with pysam.FastaFile(fasta_file) as fasta:
            # Example: Fetch reference base for a specific chromosome and position
            #chromosome = "chr17"
            #position = 4899361  # 1-based position as in VCF
            
            print("ch: " + str(chromosome) + ", " + str(position))
            # FASTA uses 0-based indexing internally
            ref_base = fasta.fetch("chr" + str(chromosome), int(position) - 10, int(position) + 10)  # Fetch one base
            
            #print(f"Reference base at {chromosome}:{position} is {ref_base}")
            return ref_base
    except Exception as e:
        return "Error: {}".format(e)

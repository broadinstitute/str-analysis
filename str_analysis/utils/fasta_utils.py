
REFERENCE_CACHE = {}

def get_reference_sequence_with_cache(fasta_obj, chrom, start_0based, end):
    """Returns the reference nucleotides in the given interval.

    Args:
        fasta_obj (object): pysam.FastaFile object
        chrom (str): chromosome name
        start_0based (int): interval start coordinate
        end (int): interval end coordinate

    Return:
        str: nucleotide sequence
    """

    if chrom not in fasta_obj:
        raise ValueError(f"Chromosome name {chrom} doesn't match any chromosome in the reference genome: {fasta_obj.keys() if hasattr(fasta_obj, 'keys') else fasta_obj.references}")

    if chrom not in REFERENCE_CACHE:
        REFERENCE_CACHE[chrom] = str(fasta_obj[chrom])

    return REFERENCE_CACHE[chrom][start_0based: end]

def get_chromosome_size_with_cache(fasta_obj, chrom):
    if chrom not in fasta_obj:
        raise ValueError(f"Chromosome name {chrom} doesn't match any chromosome in the reference genome: {fasta_obj.keys() if hasattr(fasta_obj, 'keys') else fasta_obj.references}")

    if chrom not in REFERENCE_CACHE:
        REFERENCE_CACHE[chrom] = str(fasta_obj[chrom])

    return len(REFERENCE_CACHE[chrom])


def get_reference_sequence(fasta_obj, chrom, start_1based, end_1based):
    """Returns the reference nucleotides in the given interval.

    Args:
        fasta_obj (object): pysam.FastaFile object
        chrom (str): chromosome name
        start_1based (int): interval start coordinate
        end_1based (int): interval end coordinate

    Return:
        str: nucleotide sequence
    """

    if chrom not in fasta_obj:
        raise ValueError(f"Invalid chromosome name: {chrom}")

    if fasta_obj.faidx.one_based_attributes:
        raise ValueError("The pysam.FastaFile object must be created with `one_based_attributes=False`")

    chrom_obj = fasta_obj[chrom]

    return str(chrom_obj[start_1based - 1: end_1based])


def get_chromosome_sizes(fasta_path):
    """Returns a dictionary that maps chromosome name to its size in base pairs.

    Args:
        fasta_path (str): Path of the reference genome .fasta. This function will load the corresponding .fai file to
            get the chrom. sizes.

    Return:
        dict: chromosome name to size in base pairs.
    """
    chrom_size_lookup = {}
    with open(f"{fasta_path}.fai", "rt") as fai_file:
        for line in fai_file:
            fields = line.split()
            chrom = fields[0]
            size = int(fields[1])
            chrom_size_lookup[chrom] = size
    return chrom_size_lookup


def create_normalize_chrom_function(has_chr_prefix=False):
    def normalize_chrom(chrom):
        chrom_without_prefix = chrom.replace("chr", "")
        if not has_chr_prefix:
            return chrom_without_prefix

        return f"chr{chrom_without_prefix}"

    return normalize_chrom


def create_normalize_chrom_function_for_reference_fasta(reference_fasta_path):
    import pysam
    ref_fasta = pysam.FastaFile(reference_fasta_path)
    does_chrom_start_with_chr = ref_fasta.references[0].startswith("chr")
    return create_normalize_chrom_function(does_chrom_start_with_chr)


def normalize_chrom_using_pyfaidx_fasta(pyfaidx_reference_fasta_obj, chrom):
    needs_chr_prefix = next(iter(pyfaidx_reference_fasta_obj.keys())).startswith("chr")
    chrom = chrom.replace('chr', '')
    if needs_chr_prefix:
        chrom = f"chr{chrom}"
    return chrom

def normalize_chrom_using_pysam_fasta(pysam_reference_fasta_obj, chrom):
    needs_chr_prefix = pysam_reference_fasta_obj.references[0].startswith("chr")
    chrom = chrom.replace('chr', '')
    if needs_chr_prefix:
        chrom = f"chr{chrom}"
    return chrom

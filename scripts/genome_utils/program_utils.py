"""
Utility functions for reading various file formats in SIDEARM.
This file contains functions extracted from:
- intersection_counting.py
- naive_sims.py
- interval_extraction.py
- redistribution_optimized.py
- write_encodings.py
"""

import numpy as np


def read_fasta(fasta_file):
    """
    Reads genomic sequences from a FASTA format file.
    
    Args:
        fasta_file (str): Path to the FASTA file
        
    Returns:
        dict: Dictionary mapping chromosome names to sequences
    """
    genome = {}
    with open(fasta_file, 'r') as f:
        chrom = None
        seq = ''
        for line in f:
            if line.startswith('>'):
                if chrom is not None:
                    genome[chrom] = seq
                chrom = line.strip()[1:]
                seq = ''
            else:
                seq += line.strip()
        if chrom is not None:
            genome[chrom] = seq
    return genome


def read_tf_file(tf_file):
    """
    Reads transcription factor binding sites from a file.
    
    Args:
        tf_file (str): Path to the transcription factor file
        
    Returns:
        list: List of transcription factor binding sites with parsed coordinates
    """
    tfs = []
    with open(tf_file) as f:
        for line in f:
            new_tf = line.strip().split("\t")
            new_tf[1], new_tf[2] = int(new_tf[1]), int(new_tf[2])
            tfs.append(new_tf)
    
    return tfs


def read_dmg_file(dmg_file):
    """
    Reads damage sites from a file.
    
    Args:
        dmg_file (str): Path to the damage file
        
    Returns:
        list: List of damage sites with parsed coordinates and values
    """
    dmgs = []
    with open(dmg_file) as f:
        for line in f:
            new_dmg = line.strip().split("\t")
            new_dmg[1], new_dmg[2], new_dmg[3] = int(new_dmg[1]), int(new_dmg[2]), int(new_dmg[3])
            dmgs.append(new_dmg)
    
    return dmgs


def get_intervals(interval_file):
    """
    Reads intervals from a file and assigns IDs to them.
    
    Args:
        interval_file (str): Path to the interval file
        
    Returns:
        list: List of intervals with parsed coordinates and assigned IDs
    """
    intervalId = 0
    intervals = []
    with open(interval_file) as f:
        for line in f:
            interval = line.strip().split("\t")
            interval.append(intervalId)
            interval[1], interval[2] = int(interval[1]), int(interval[2])
            intervals.append(interval)
            intervalId += 1
    
    return intervals


def hold_interval_indices(interval_file):
    """
    Reads the start positions of intervals from a file.
    
    Args:
        interval_file (str): Path to the interval file
        
    Returns:
        list: List of tuples (chromosome, start_position)
    """
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            
            start = int(cur[1])
            intervals.append((cur[0], start))
    
    return intervals


def collapse_dmg_sites(input_file):
    """
    Collapses damage sites from a BED file format.
    
    Args:
        input_file (str): Path to the input damage file
        
    Returns:
        list: List of consolidated damage sites
    """
    # Initialize variables to keep track of the current damage site
    current_chrom = None
    current_start = None
    current_end = None
    current_strand = None
    count = 0
    matched = 0

    closed = True

    # Initialize an empty list to store the consolidated sites
    consolidated_sites = []

    with open(input_file, 'r') as infile:
        for line in infile:
            chrom, start, end, strand = line.strip().split()
            start, end = int(start), int(end)

            if strand == "+":
                strand = "-"
            elif strand == "-":
                strand = "+"

            # Check if this is the first line
            if closed == True:
                current_chrom, current_start, current_end, current_strand = chrom, start, end, strand
                count = 1
                matched = 0
                closed = False

            else:
                # If the current line is adjacent to the previous and in the same strand
                if chrom == current_chrom and strand == current_strand and start == current_start:
                    count += 1

                else:
                    # Save the previous range to the consolidated_sites list
                    matched += 1
                    if count == matched:
                        consolidated_sites.append((current_chrom, current_start, current_end, count, current_strand))
                        closed = True

        # Don't forget to save the last range after finishing the loop
        consolidated_sites.append([current_chrom, current_start, current_end, count, current_strand])

    return consolidated_sites


def rev_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    
    Args:
        seq (str): DNA sequence
        
    Returns:
        str: Reverse complement of the input sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def decode(site_code, interval_map):
    """
    Decode a site code into genomic coordinates.
    
    Args:
        site_code (int): Encoded site code
        interval_map (list): Mapping of intervals
        
    Returns:
        tuple: Decoded genomic coordinates
    """
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', f'{int(string[-6:-4], 2)}{int(string[-4:-2], 2)}{int(string[-2:], 2)}')

    string = string[1:]

    region = int(string[:17], 2)
    pos = int(string[17:26], 2)

    chrom, start = interval_map[region]
    pos += start

    strand = int(string[26])
    match strand:
        case 0:
            strand = '+'
        case 1:
            strand = '-'
    
    dmg = int(string[27:31], 2)

    return (chrom, pos, pos+1, dmg, strand)


def decompress_mut(site_code, interval_map):
    """
    Decompress a mutation site code into genomic coordinates.
    
    Args:
        site_code (int): Encoded mutation site code
        interval_map (list): Mapping of intervals
        
    Returns:
        tuple: Decompressed mutation information
    """
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', int(string[-4:], 2))

    string = string[1:]

    region = int(string[:16], 2)
    pos = int(string[16:25], 2)
    dmg = int(string[25:], 2)

    chrom, start = interval_map[region]
    pos += start

    return (chrom, pos, dmg)


def decode_nochrom(site_code):
    """
    Decode a site code without chromosome information.
    
    Args:
        site_code (int): Encoded site code
        
    Returns:
        tuple: Decoded site information
    """
    binary_str = f'{site_code:032b}'

    # If the first bit is 1 => "NEW"
    if int(binary_str[0]) == 1:
        return ('NEW', int(binary_str[-6:], 2))

    # Otherwise, skip the first 2 bits
    binary_str = binary_str[2:]

    region = int(binary_str[:17], 2)
    pos = int(binary_str[17:26], 2)
    strand_bit = int(binary_str[26])
    strand = '+' if strand_bit == 0 else '-'
    dmg = int(binary_str[27:30], 2)

    return (region, pos, strand, dmg)


def pack_tuple_nochrom(tup):
    """
    Re-encode (region, pos, strand, damage) into a 32-bit int.
    
    Args:
        tup (tuple): Site information to encode
        
    Returns:
        int: Packed site code
    """
    region_val, pos_val, strand_sym, dmg_val = tup

    # Basic bounds checks
    if not (0 <= region_val < 2**17):
        raise ValueError(f"[pack_tuple_nochrom] Region {region_val} >= 2^17?")
    if not (0 <= pos_val < 2**9):
        raise ValueError(f"[pack_tuple_nochrom] Position {pos_val} >= 2^9?")
    if not (0 <= dmg_val < 2**4):
        raise ValueError(f"[pack_tuple_nochrom] Damage {dmg_val} >= 2^4?")

    region_bits = f"{region_val:017b}"
    pos_bits    = f"{pos_val:09b}"
    strand_bit  = '0' if strand_sym == '+' else '1'
    dmg_bits    = f"{dmg_val:04b}"

    # Combine the bits
    bit_str = region_bits + pos_bits + strand_bit + dmg_bits
    int_val = int(bit_str, 2)
    return int_val


def decompress_3mer(site_code, region_bits=16, pos_bits=9, strand_bits=0):
    """
    Decompress a 3-mer site code.
    
    Args:
        site_code (int): Encoded 3-mer site code
        region_bits (int): Number of bits for region
        pos_bits (int): Number of bits for position
        strand_bits (int): Number of bits for strand
        
    Returns:
        tuple: Decompressed site information
    """
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', int(string[-4:], 2))

    string = string[1:]

    dmg_bits = 31 - region_bits - pos_bits - strand_bits
    region = int(string[:region_bits], 2)
    pos = int(string[region_bits:region_bits+pos_bits], 2)
    
    if strand_bits > 0:
        strand = int(string[region_bits+pos_bits:region_bits+pos_bits+strand_bits], 2)
        dmg = int(string[region_bits+pos_bits+strand_bits:], 2)
        return (region, pos, strand, dmg)
    else:
        dmg = int(string[region_bits+pos_bits:], 2)
        return (region, pos, dmg)


def compress_3mer(tup, region_bits=16, pos_bits=9, strand_bits=0):
    """
    Compress a 3-mer site tuple into an encoded value.
    
    Args:
        tup (tuple): Site information to encode
        region_bits (int): Number of bits for region
        pos_bits (int): Number of bits for position
        strand_bits (int): Number of bits for strand
        
    Returns:
        int: Compressed site code
    """
    dmg_bits = 31 - region_bits - pos_bits - strand_bits
    region = f'{tup[0]:0{region_bits}b}'
    pos = f'{tup[1]:0{pos_bits}b}'
    if len(tup) > 3:
        strand = f'{tup[2]:0{strand_bits}b}'
        dmg = f'{tup[3]:0{dmg_bits}b}'
    else:
        strand = ''

        if tup[2] >= 2**dmg_bits:
            print(f'{tup[2]} is too large for {dmg_bits} bits')
            dmg = f'{2**dmg_bits - 1:0{dmg_bits}b}'
        else:
            dmg = f'{tup[2]:0{dmg_bits}b}'

    return int(region + pos + strand + dmg, 2)


def match_cpd_4mer(seq):
    """
    Create a numerical representation for each 4mer that can be easily interpreted.
    
    Args:
        seq (str): 4-mer DNA sequence
        
    Returns:
        tuple: Numerical representation and strand
    """
    seq = seq.upper()

    if len(seq) != 4:
        print("broken")
        return 64, "."

    if 'T' in seq[1:3] or 'C' in seq[1:3]:
        strand = "+"
    
    else:
        strand = "-"
        seq = rev_complement(seq)

    f1 = seq[0]
    dimer = seq[1:3]
    f2 = seq[3]

    tri = ""

    match f1:
        case "A":
            tri += "0"
        
        case "C":
            tri += "1"

        case "G":
            tri += "2"
        
        case "T":
            tri += "3"
        
        case _:
            tri += "1000"

    match dimer:
        case "CC":
            tri += "0"
        
        case "CT":
            tri += "1"

        case "TC":
            tri += "2"
        
        case "TT":
            tri += "3"

        case _:
            tri += "1000"
    
    match f2:
        case "A":
            tri += "0"
        
        case "C":
            tri += "1"

        case "G":
            tri += "2"
        
        case "T":
            tri += "3"
        
        case _:
            tri += "1000"
    
    return int(tri, 4), strand


def match_bpde_3mer(seq):
    """
    Create a numerical representation for BPDE 3-mers.
    
    Args:
        seq (str): 3-mer DNA sequence
        
    Returns:
        tuple: Numerical representation and strand
    """
    seq = seq.upper()

    if len(seq) != 3 or seq[1] not in ['G', 'C']:
        return 16, "."

    if seq[1] == 'G':
        strand = "+"
    
    else:
        strand = "-"
        seq = rev_complement(seq)
    

    tri = ""
    match seq[0]:
        case "A":
            tri += "0"
        
        case "C":
            tri += "1"

        case "G":
            tri += "2"
        
        case "T":
            tri += "3"
        
        case _:
            tri += "100"
    
    match seq[2]:
        case "A":
            tri += "0"
        
        case "C":
            tri += "1"

        case "G":
            tri += "2"
        
        case "T":
            tri += "3"

        case _:
            tri += "100"
    
    return int(tri, 4), strand


def match_C_3mer(seq):
    """
    Match C-centered 3-mers and encode them.
    
    Args:
        seq (str): 3-mer DNA sequence
        
    Returns:
        int: Encoded value of the 3-mer
    """
    mapping_edge = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    seq = seq.upper()

    if any(base not in 'ACGT' for base in seq) or len(seq) != 3 or seq[1] not in ['C', 'G']:
        return 32

    if seq[1] != 'C':
        seq = rev_complement(seq)

    if seq[1] != 'C':
        return 32
        
    try:
        return int(f'{mapping_edge[seq[0]]}{mapping_edge[seq[2]]}', 2)
    except KeyError:
        return 32


def redistribute(lists, total_items):
    """
    Redistribute the total damage within a sequence to its components.
    
    Args:
        lists (list): List of positions and values
        total_items (int): Total number of items to distribute
        
    Returns:
        list: Redistributed damage values
    """
    try:
        total_items = int(total_items)  # Ensure total_items is an integer
    except (ValueError, TypeError):
        raise ValueError(f"total_items must be convertible to integer, got: {total_items}")

    if not lists:
        raise ValueError("Input lists cannot be empty")

    # Unpack the lists with type checking
    chrom = [row[0] for row in lists]
    pos_1 = [int(row[1]) for row in lists]  # Ensure positions are integers
    pos_2 = [int(row[2]) for row in lists]
    dmg = [float(row[3]) for row in lists]
    strand = [row[4] for row in lists]
    
    intersect = np.array([pos_1, dmg]).T

    num_destinations = len(intersect[:, 0])
    if num_destinations == 0:
        raise ValueError("No valid destinations found")

    # Define the probabilities of selecting each destination
    probabilities = np.ones(num_destinations) / num_destinations

    # Distribute the items randomly
    destinations = np.random.choice(num_destinations, total_items, p=probabilities)

    # Count the number of items per destination
    counts = np.bincount(destinations, minlength=num_destinations)

    # Convert counts to integer list
    counts = counts.astype(int).tolist()

    # Return the combined lists
    return [list(x) for x in zip(chrom, pos_1, pos_2, counts, strand)]
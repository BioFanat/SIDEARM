import pandas as pd
from intervaltree import IntervalTree, Interval
import struct
import numpy as np
import math

def read_fasta(fasta_file):
    """Reads genomic sequences in .fasta format
    
    Args:
        fasta_file: Path to fasta file
        
    Returns:
        dict: Dictionary of chromosome sequences
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
        genome[chrom] = seq
    return genome

def rev_complement(seq):
    """Return the reverse complement of a DNA sequence
    
    Args:
        seq: DNA sequence
        
    Returns:
        str: Reverse complement of input sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def comp_strands(dmg_strand, tf_strand):
    """Compare damage and TF strands to see if damage was on GGAA (True) or TTCC (False)
    
    Args:
        dmg_strand: Damage strand
        tf_strand: TF strand
        
    Returns:
        bool: True if strands are the same, indicating "GGAA", False otherwise
    """
    return dmg_strand == tf_strand  # if the strands are the same, then this is "GGAA" otherwise "TTCC"

def match_cpd_4mer(seq):
    """Create a numerical representation for each 4mer that can be easily interpreted
    
    Args:
        seq: 4mer DNA sequence
        
    Returns:
        tuple: (index, strand) - Index is an integer representation, strand is +/-
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
    """Create a numerical representation for each 3mer that can be easily interpreted
    
    Args:
        seq: 3mer DNA sequence
        
    Returns:
        tuple: (index, strand) - Index is an integer representation, strand is +/-
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

def read_tf_file(tf_file):
    """Read a TF binding site file
    
    Args:
        tf_file: Path to TF binding site file
        
    Returns:
        list: List of TF binding sites
    """
    tfs = []
    with open(tf_file) as f:
        for line in f:
            new_tf = line.strip().split("\t")
            new_tf[1], new_tf[2] = int(new_tf[1]), int(new_tf[2])
            tfs.append(new_tf)
    
    return tfs

def read_dmg_file(dmg_file):
    """Read a damage file
    
    Args:
        dmg_file: Path to damage file
        
    Returns:
        list: List of damage sites
    """
    dmgs = []
    with open(dmg_file) as f:
        for line in f:
            new_dmg = line.strip().split("\t")
            new_dmg[1], new_dmg[2], new_dmg[3] = int(new_dmg[1]), int(new_dmg[2]), int(new_dmg[3])
            dmgs.append(new_dmg)
    
    return dmgs

def hold_interval_indices(interval_file):
    """Hold interval indices for efficient lookup
    
    Args:
        interval_file: Path to interval file
        
    Returns:
        list: List of (chrom, start) tuples
    """
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            
            start = int(cur[1])
            intervals.append((cur[0], start))
    
    return intervals

def process_intervals(tfbs, radius):
    """Process TFBS intervals and create an interval tree
    
    Args:
        tfbs: List of TFBS
        radius: Radius around midpoint to create intervals
        
    Returns:
        dict: Dictionary of IntervalTree objects by chromosome
    """
    tf_intervals = {}

    # Global motif_length - needs to be set or passed in real application
    motif_length = tfbs[0][2] - tfbs[0][1]

    for tf in tfbs:
        chrom = tf[0]
        midpoint = (tf[1] + tf[2] - 1) / 2
        tf_strand = tf[3]

        match tf_strand:
            case "+":
                midpoint = int(math.floor(midpoint))
            
            case "-":
                midpoint = int(math.ceil(midpoint))

        if chrom not in tf_intervals:  # tf[0] is the chromosome
            tf_intervals[chrom] = IntervalTree()
            tf_intervals[chrom].add(Interval(midpoint - radius, midpoint + radius, tf_strand))
        else:
            tf_intervals[chrom].add(Interval(midpoint - radius, midpoint + radius, tf_strand))
    
    return tf_intervals

def overlap(damage, intervals, genome):
    """Check if a damage site overlaps with any TFBS
    
    Args:
        damage: Damage site (chrom, start, end, value, strand)
        intervals: Dictionary of interval trees by chromosome
        genome: Genome sequence dictionary
        
    Returns:
        list: List of overlaps (chrom, distance, damage_value, strand_comparison)
    """
    chrom = damage[0]
    pos = (damage[1], damage[2])  # For CPD
    dmg_val = damage[3]
    dmg_strand = damage[4]

    seq = None
    
    if chrom in genome:
        seq = genome[chrom][pos[0]-1:pos[1]+2]  # For CPD

    code, strand = match_cpd_4mer(seq)  # FOR CPD
    if not seq or code > 64:
        return []
    
    overlaps = []
    if chrom in intervals:
        overlapping = intervals[chrom].overlap(*pos)  # FOR CPD
        for overlap in overlapping:
            mid = (overlap[0] + overlap[1]) // 2
            tf_strand = overlap[2]
            
            dmg_marker = damage[2] if tf_strand == "-" else damage[1]  # FOR CPD
            dist = int(dmg_marker - mid)
            
            # Optional: get sequence context
            sequence = genome[chrom][mid-15:mid+16].upper()

            if tf_strand == "-":
                dist = -dist

            overlaps.append((chrom, dist, dmg_val, comp_strands(dmg_strand, tf_strand)))
    
    return overlaps

def overlap_muts(mutation, intervals, genome):
    """Check if a mutation site overlaps with any TFBS
    
    Args:
        mutation: Mutation site (chrom, pos)
        intervals: Dictionary of interval trees by chromosome
        genome: Genome sequence dictionary
        
    Returns:
        list: List of overlaps (chrom, distance)
    """
    chrom, pos = mutation[0], mutation[1]

    overlaps = []
    if chrom in intervals:
        overlapping = intervals[chrom].at(pos)
        for overlap in overlapping:
            mid = (overlap[0] + overlap[1]) // 2
            tf_strand = overlap[2]
            
            dist = int(pos - mid)

            # Optional: get sequence context
            sequence = genome[chrom][mid-15:mid+16].upper()

            if tf_strand == "-":
                dist = -dist

            overlaps.append((chrom, dist))
    
    return overlaps

def aggregate_overlaps(overlaps, window, num_tfs):
    """Aggregate overlap counts by position and strand
    
    Args:
        overlaps: List of overlaps
        window: Window size around midpoint
        num_tfs: Number of TFs
        
    Returns:
        DataFrame: Aggregated overlap counts by position and strand
    """
    rel_dmg_same = []
    rel_dmg_diff = []
    for i in range(2*window + 1):
        rel_dmg_same.append([i - window + 0.5, 0.0, "Same"])
        rel_dmg_diff.append([i - window + 0.5, 0.0, "Diff"])
    
    for overlap in overlaps:
        dist = overlap[1]
        if abs(dist) <= window:    
            dmg = overlap[2]
            if overlap[3]:
                rel_dmg_same[dist + window][1] += dmg
            else:
                rel_dmg_diff[dist + window][1] += dmg

    rel_dmg = rel_dmg_same + rel_dmg_diff
    rel = pd.DataFrame(rel_dmg, columns = ["Pos", "Dmg", "Strand"])
    rel = rel.set_index("Pos")
    
    return rel

def process_dmgs(plus_dmg_file, minus_dmg_file, tree, window, genome):
    """Process damages from plus and minus strands
    
    Args:
        plus_dmg_file: Path to plus strand damage file
        minus_dmg_file: Path to minus strand damage file
        tree: Dictionary of interval trees by chromosome
        window: Window size around midpoint
        genome: Genome sequence dictionary
        
    Returns:
        DataFrame: Processed damages with position and strand information
    """
    rel_dmg_same = []
    rel_dmg_diff = []
    for i in range(2*window + 1):
        rel_dmg_same.append([i - window + 0.5, 0.0, "Same"])
        rel_dmg_diff.append([i - window + 0.5, 0.0, "Diff"])

    with open(plus_dmg_file) as f:
        for line in f:
            new_dmg = line.strip().split("\t")
            new_dmg[1], new_dmg[2], new_dmg[3] = int(new_dmg[1]), int(new_dmg[2]), int(new_dmg[3])

            if new_dmg[3] > 0:
                for match in overlap(new_dmg, tree, genome):
                    dist = match[1]
                    if abs(dist) <= window:    
                        dmg = match[2]
                        if match[3]:
                            rel_dmg_same[dist + window][1] += dmg
                        else:
                            rel_dmg_diff[dist + window][1] += dmg
    
    with open(minus_dmg_file) as f:
        for line in f:
            new_dmg = line.strip().split("\t")
            new_dmg[1], new_dmg[2], new_dmg[3] = int(new_dmg[1]), int(new_dmg[2]), int(new_dmg[3])
            if new_dmg[3] > 0:
                for match in overlap(new_dmg, tree, genome):
                    dist = match[1]
                    if abs(dist) <= window:    
                        dmg = match[2]
                        if match[3]:
                            rel_dmg_same[dist + window][1] += dmg
                        else:
                            rel_dmg_diff[dist + window][1] += dmg
    rel_dmg = rel_dmg_same + rel_dmg_diff
    rel = pd.DataFrame(rel_dmg, columns = ["Pos", "Dmg", "Strand"])
    rel = rel.set_index("Pos")

    return rel

def process_muts(mutation_file, tree, window, genome):
    """Process mutations
    
    Args:
        mutation_file: Path to mutation file
        tree: Dictionary of interval trees by chromosome
        window: Window size around midpoint
        genome: Genome sequence dictionary
        
    Returns:
        DataFrame: Processed mutations with position information
    """
    rel_mut = []

    for i in range(1, 2*window):
        rel_mut.append([i - window, 0.0])

    with open(mutation_file) as f:
        for line in f:
            new_mut = line.strip().split("\t")
            new_mut[1], new_mut[2] = int(new_mut[1]), int(new_mut[2])

            for match in overlap_muts(new_mut, tree, genome):
                dist = match[1]
                if abs(dist) < window:
                    rel_mut[dist + window - 1][1] += 1
    
    rel = pd.DataFrame(rel_mut, columns = ["Pos", "Mut"])
    rel = rel.set_index("Pos")

    return rel

def redistribute(lists, total_items):
    """Randomly redistribute items across positions
    
    Args:
        lists: List of sites
        total_items: Total number of items to redistribute
        
    Returns:
        list: List with redistributed items
    """
    # Type checking and conversion
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

def decode(site_code, interval_map):
    """Decode a 32-bit integer site code to genomic coordinates
    
    Args:
        site_code: 32-bit integer site code
        interval_map: List of (chrom, start) tuples
        
    Returns:
        tuple: Decoded site information
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
    """Decompress mutation site code
    
    Args:
        site_code: 32-bit integer site code
        interval_map: List of (chrom, start) tuples
        
    Returns:
        tuple: Decompressed mutation site information
    """
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', int(string[-5:], 2))

    string = string[1:]

    region = int(string[:17], 2)
    pos = int(string[17:26], 2)
    dmg = int(string[26:], 2)

    chrom, start = interval_map[region]
    pos += start

    return (chrom, pos, dmg)

def pack_tuple_nochrom(tup):
    """Pack a tuple without chromosome information into a 32-bit integer
    
    Args:
        tup: Tuple of (region, position, strand, damage)
        
    Returns:
        int: Packed 32-bit integer
    """
    if not (0 <= tup[0] < 2**17):  # Check region fits in 17 bits
        raise ValueError(f"Region {tup[0]} too large for 17 bits")
    if not (0 <= tup[1] < 2**9):   # Check position fits in 9 bits
        raise ValueError(f"Position {tup[1]} too large for 9 bits")
    if not (0 <= tup[3] < 2**3):   # Check damage fits in 3 bits
        raise ValueError(f"Damage {tup[3]} too large for 3 bits")
        
    region = f'{tup[0]:017b}'
    pos = f'{tup[1]:09b}'
    strand = '0' if tup[2] == '+' else '1'
    dmg = f'{tup[3]:03b}'
    
    packed = int(region + pos + strand + dmg, 2)
    assert packed < 2**32, f"Packed value {packed} exceeds 32 bits"
    return packed
import pandas as pd
from intervaltree import IntervalTree, Interval
import struct
import time
from multiprocessing import Pool
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
import logomaker as lm
import statsmodels.stats.multitest as smstats
import scipy.stats as stats
import argparse
import math
import psutil
import os

def get_memory_usage(location=""):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 * 1024)  # Convert to MB
    print(f"[MEMORY CHECK - {location}] Process {os.getpid()}: {mem:.2f} MB")

def read_fasta(fasta_file): # reads genomic sequences - in .fasta format

    get_memory_usage("Before genome load")

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

    get_memory_usage("After genome load")

    return genome

def process_tfbs_file(input_file, output_file): # reads TF damages - in .bed format
    out = open(output_file, 'w')
    with open(input_file) as f:
        for line in f:
            fields = line.strip().split(" ")
            chrom = fields[0]
            start = int(fields[1])-1
            end = int(fields[2])-1
            strand = fields[5]
            out.write(f'{chrom}\t{start}\t{end}\t{strand}\n')
        
        out.close()

def read_tf_file(tf_file):
    tfs = []
    with open(tf_file) as f:
        for line in f:
            new_tf = line.strip().split("\t")
            new_tf[1], new_tf[2] = int(new_tf[1]), int(new_tf[2])
            tfs.append(new_tf)
    
    return tfs
    

def read_dmg_file(dmg_file):
    dmgs = []
    with open(dmg_file) as f:
        for line in f:
            new_dmg = line.strip().split("\t")
            new_dmg[1], new_dmg[2], new_dmg[3] = int(new_dmg[1]), int(new_dmg[2]), int(new_dmg[3])
            # print(new_dmg)
            dmgs.append(new_dmg)
    
    return dmgs

def get_intervals(interval_file):
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

def extend_intervals(interval_file, output_file, radius):
    out = open(output_file, 'w')
    with open(interval_file) as f:
        for line in f:
            interval = line.strip().split("\t")
            interval[1], interval[2] = str(int(interval[1])-radius), str(int(interval[2])+radius)
            out.write("\t".join(interval) + "\n")
    
        out.close()

def match_intervals(tfs, intervals, output_file):
    new_tfs = []
    tf_point = 0
    interval_point = 0

    out = open(output_file, 'w')

    while tf_point < len(tfs) and interval_point < len(intervals):
        cur_interval = intervals[interval_point]
        cur_chr = cur_interval[0]
        cur_int = (cur_interval[1], cur_interval[2])
        cur_id = cur_interval[3]

        while tf_point < len(tfs) and tfs[tf_point][0] == cur_chr and (tfs[tf_point][1] >= cur_int[0] and tfs[tf_point][2] <= cur_int[1]):
            cur_tf = tfs[tf_point]
            cur_tf.append(cur_id)

            out.write("\t".join([str(c) for c in cur_tf]) + "\n")
            new_tfs.append(cur_tf)
            tf_point += 1

        interval_point += 1
    
    # print(len(tfs), len(intervals), len(new_tfs), tf_point, interval_point)
    out.close()

def comp_strands(dmg_strand, tf_strand): # compare the damage and tf strands to see if the damage was on GGAA (true) or TTCC (false)
    return dmg_strand == tf_strand # if the strands are the same, then this is "GGAA" otherwise "TTCC"

def decode_ind(index):
    bin_rep = f'{index:06b}'
    first = bin_rep[:2]
    middle = bin_rep[2:4]
    last = bin_rep[4:]

    nuc = ['A', 'C', 'G', 'T']
    dimer = ['CC', 'CT', 'TC', 'TT']

    return f'{nuc[int(first, 2)]}{dimer[int(middle, 2)]}{nuc[int(last, 2)]}'

def match_cpd_4mer(seq): #create a numerical representation for each 4mer that can be easily interpreted

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

def rev_complement(seq): # find the reverse complement for a sequence
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def process_intervals(tfbs, radius):

    get_memory_usage("Starting interval tree creation")

    tf_intervals = {}

    global motif_length
    motif_length = tfbs[0][2]-tfbs[0][1]
    #print(motif_length, tfbs[0])

    for tf in tfbs:
        chrom = tf[0]
        midpoint = (tf[1] + tf[2] - 1) / 2
        tf_strand = tf[3]

        match tf_strand:
            case "+":
                midpoint = int(math.floor(midpoint))
            
            case "-":
                midpoint = int(math.ceil(midpoint))

        if chrom not in tf_intervals: #tf[0] is the chromosome
            tf_intervals[chrom] = IntervalTree()
            tf_intervals[chrom].add(Interval(midpoint - radius, midpoint + radius, tf_strand))
        else:
            tf_intervals[chrom].add(Interval(midpoint - radius, midpoint + radius, tf_strand))
    

    get_memory_usage("After interval tree creation")
    return tf_intervals

def aggregate_overlaps(overlaps, window, num_tfs):
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
    #rel["Dmg"] = rel["Dmg"].div(num_tfs)
    return rel

def intersection(tf_file, pd_file, md_file):

    seen = set()

    with open(tf_file) as f:
        curInd = 0

        for line in f:
            new_tf = line.strip().split("\t")

    
    pd = open(pd_file, 'r')
    md = open(md_file, 'r')

def hold_interval_indices(interval_file):
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            
            start = int(cur[1])
            intervals.append((cur[0], start))
    
    return intervals
    
def decode(site_code, interval_map):

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
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', int(string[-4:], 2))

    string = string[1:]

    region = int(string[:16], 2)
    pos = int(string[16:25], 2)
    dmg = int(string[25:], 2)

    # if region >= len(interval_map):
    #     print(region, pos)
    #     chrom, start = 'chr1 ', 0
   
    chrom, start = interval_map[region]
    pos += start

    return (chrom, pos, dmg)

def overlap(damage, intervals):
    chrom = damage[0]
    pos = (damage[1], damage[2]) # - FOR CPD
    #pos = damage[1] # FOR BPDE
    dmg_val = damage[3]
    dmg_strand = damage[4]

    seq = None
    
    #global genome
    if chrom in genome:
        seq = genome[chrom][pos[0]-1:pos[1]+2] #- FOR CPD
        #seq = genome[chrom][pos-1:pos+2] # FOR BPDE

    code, strand = match_cpd_4mer(seq) # FOR CPD
    #code, strand = match_bpde_3mer(seq) # FOR BPDE
    if not seq or code > 64: #or code[1] == '3':
        return []
    
    # dmg_marker = damage[1]
    # if dmg_strand == "+":
    #     dmg_marker = damage[2]

    overlaps = []
    if chrom in intervals:
        overlapping = intervals[chrom].overlap(*pos) #- FOR CPD
        #overlapping = intervals[chrom].at(pos) # FOR BPDE
        for overlap in overlapping:
            mid = (overlap[0] + overlap[1]) // 2
            tf_strand = overlap[2]
            
            dmg_marker = damage[2] if tf_strand == "-" else damage[1]# - FOR CPD
            #dmg_marker = pos
            dist = int(dmg_marker - mid)
            

            #global genome
            sequence = genome[chrom][mid-15:mid+16].upper()

            if tf_strand == "-":
                dist = -dist

                global minus_tf
                minus_tf.append(sequence)

            else:

                global plus_tf
                plus_tf.append(sequence)

            overlaps.append((chrom, dist, dmg_val, comp_strands(dmg_strand, tf_strand)))
    
    return overlaps

def overlap_muts(mutation, intervals):
    chrom, pos = mutation[0], mutation[1]
    #pos = (mutation[1], mutation[2])

    #print(mutation)

    # dmg_marker = damage[1]
    # if dmg_strand == "+":
    #     dmg_marker = damage[2]

    overlaps = []
    if chrom in intervals:
        overlapping = intervals[chrom].at(pos)
        for overlap in overlapping:
            #print(f'{chrom}\t{pos}\t{overlap}')
            mid = (overlap[0] + overlap[1]) // 2
            tf_strand = overlap[2]
            
            #dmg_marker = damage[2] if tf_strand == "-" else damage[1]
            dist = int(pos - mid)

            #global genome
            sequence = genome[chrom][mid-15:mid+16].upper()

            if tf_strand == "-":
                dist = -dist

                global minus_tf
                minus_tf.append(sequence)

            else:

                global plus_tf
                plus_tf.append(sequence)

            overlaps.append((chrom, dist))
    
    return overlaps

def process_run(file_path, intervals, tree, window):
    rel_dmg_same = []
    rel_dmg_diff = []
    for i in range(2*window + 1):
        rel_dmg_same.append(0.0)
        rel_dmg_diff.append(0.0)
    
    bin_file = open(file_path, 'rb')
    data = bin_file.read(4)

    # cur_4mer = '000'
    first_run = True

    while data:
        # next_data = bin_file.read(4)
        source = decode(struct.unpack('I', data)[0], intervals)

        if source[0] == "NEW":
            #cur_4mer += 1
            # print("NEW 4MER", source[1])
            cur_4mer = source[1]
            # if cur_4mer[1] == '3':
            #     print('non_mutagenic')

            if source[1] == '000' and not first_run:
                break
            
            first_run = False

        else:
            if source[3] > 0:    
                for match in overlap(source, tree):
                    dist = match[1]
                    # if dist == 1:
                    #     print(decode_ind(cur_4mer), source)

                    if abs(dist) <= window:    
                        dmg = match[2]
                        if match[3]:
                            rel_dmg_same[dist + window] += dmg
                        else:
                            rel_dmg_diff[dist + window] += dmg
        
        data = bin_file.read(4)
                
    bin_file.close()

    return rel_dmg_same + rel_dmg_diff

def process_dmgs(plus_dmg_file, minus_dmg_file, tree, window):
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
                for match in overlap(new_dmg, tree):
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
                for match in overlap(new_dmg, tree):
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

def process_muts(mutation_file, tree, window):
    rel_mut = []

    for i in range(1, 2*window):
        rel_mut.append([i - window, 0.0])

    with open(mutation_file) as f:
        for line in f:
            new_mut = line.strip().split("\t")
            new_mut[1], new_mut[2] = int(new_mut[1]), int(new_mut[2])

            for match in overlap_muts(new_mut, tree):
                dist = match[1]
                if abs(dist) < window:
                    rel_mut[dist + window - 1][1] += 1
    
    rel = pd.DataFrame(rel_mut, columns = ["Pos", "Mut"])
    rel = rel.set_index("Pos")

    return rel

def process_mut_run(file_path, intervals, tree, window):
    rel_mut = []

    for i in range(1, 2*window):
        rel_mut.append(0.0)
    
    bin_file = open(file_path, 'rb')
    data = bin_file.read(4)

    #cur_4mer = 0

    while data:
        data = bin_file.read(4)
        source = decompress_mut(struct.unpack('I', data)[0], intervals)

        if source[0] == "NEW":
            if source[1] == 0:
                break

        else:
            if source[2] > 0:    
                # if source[2] > 1:
                #     print("glitch")
                for match in overlap_muts(source, tree):
                    dist = match[1]
                    if abs(dist) < window:
                        rel_mut[dist + window - 1] += source[2]
        
                
    bin_file.close()

    return rel_mut

def perform_run(intervals, tree, radius, file_path_id):
    file_path = f'/usr/xtmp/bc301/sim_atac_nb/acc_run_{file_path_id}.bin'
    return process_run(file_path, intervals, tree, radius)

def build_runs(num_runs, range_map, tf_tree, window):
    total_runs = []
    func = partial(perform_run, range_map, tf_tree, window)

    with Pool() as p:
        total_runs = p.map(func, range(1, num_runs+1))

    # for i in range(1, num_runs+1):
    #     total_runs.append(func(i))
    
    df = pd.DataFrame(total_runs, columns = [i - window + 0.5 for i in range(2*window + 1)] * 2)
    df = df.transpose()

    df.index.name = "Pos"
    df.columns = [f'Run {i}' for i in range(1, num_runs+1)]
    df[['min', 'max', 'median', 'mean']] = df.apply(lambda row: pd.Series({'min': np.min(row), 'max': np.max(row), 'median': np.median(row), 'mean': np.mean(row)}), axis=1)
    
    return df

def perform_mut(intervals, tree, radius, mode, file_path_id):

    process = psutil.Process(os.getpid())
    threads = process.threads()
    
    print(f"\n=== Process {os.getpid()} Information ===")
    print(f"Number of threads: {len(threads)}")
    print("Thread IDs:", [thread.id for thread in threads])
    
    # Get memory info for this process
    mem = process.memory_info().rss / (1024 * 1024)
    print(f"Process memory: {mem:.2f} MB")

    paths = {'actual': f'/usr/xtmp/bc301/sim_data_skin_mut_atac/acc_run_{file_path_id}.bin', 'potential': f'/usr/xtmp/bc301/sim_data_potential_mut/acc_run_{file_path_id}.bin'}
    # file_path = f'/usr/xtmp/bc301/sim_data_mut_cor/acc_run_{file_path_id}.bin'
    # file_path = f'/usr/xtmp/bc301/sim_data_potential_mut/acc_run_{file_path_id}.bin'

    mem_after = process.memory_info().rss / (1024 * 1024)
    print(f"Process {os.getpid()} final memory: {mem_after:.2f} MB")
    
    return process_mut_run(paths[mode], intervals, tree, radius)

def build_muts(num_runs, range_map, tf_tree, window, mode):
    total_runs = []
    func = partial(perform_mut, range_map, tf_tree, window, mode)

    with Pool() as p:
        total_runs = p.map(func, range(1, num_runs+1))
    
    df = pd.DataFrame(total_runs, columns = [i - window for i in range(1, 2*window)])
    df = df.transpose()

    df.index.name = "Pos"
    df.columns = [f'Run {i}' for i in range(1, num_runs+1)]
    df[['min', 'max', 'median', 'mean']] = df.apply(lambda row: pd.Series({'min': np.min(row), 'max': np.max(row), 'median': np.median(row), 'mean': np.mean(row)}), axis=1)
    
    return df

def add_stats(data, num_runs):
    data.reset_index(inplace=True) 
    
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)

    # Create mask for flanks
    flank_mask = ((data['Pos'] < motif_range[0] - 5) | (data['Pos'] > motif_range[1] + 5))
    
    # Get run columns with correct slice
    run_cols = slice(3, -4)  # Changed from -3 to -4 to account for all summary statistics
    
    if flank_mask.any():
        # Calculate factor using masked data
        factor = (data.loc[flank_mask, 'Dmg'].mean() / 
                 data.loc[flank_mask, 'mean'].mean())
        print(f'Scaling Factor: {factor}')
        
        # Modify run columns in place
        data.iloc[:, run_cols] *= factor
    
    # Get runs view after modification
    runs = data.iloc[:, run_cols]
    exp = data['Dmg']
        
    # Calculate p-values
    bool_comparison_top = runs.gt(exp, axis=0)
    bool_comparison_bottom = runs.lt(exp, axis=0)

    # Calculate FDR-corrected p-values
    data['P-Value-Top'] = smstats.fdrcorrection(
        bool_comparison_top.sum(axis=1) / num_runs)[1]
    data['P-Value-Bottom'] = smstats.fdrcorrection(
        bool_comparison_bottom.sum(axis=1) / num_runs)[1]

    # Recalculate all summary statistics after scaling
    data['min'] = runs.min(axis=1)
    data['max'] = runs.max(axis=1)
    data['median'] = runs.median(axis=1)
    data['mean'] = runs.mean(axis=1)
    data['std'] = runs.std(axis=1)
    data['z'] = (data['Dmg'] - data['mean']) / (data['std'] + 1e-9)

def add_stats_split(data, num_runs):
    data.reset_index(inplace=True)
    
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)
    
    # Create masks for flanks (used for calculating scaling factors)
    left_flank_mask = (data['Pos'] < motif_range[0] - 5)
    right_flank_mask = (data['Pos'] > motif_range[1] + 5)
    
    # Create masks for applying scaling (entire left/right halves)
    left_half_mask = (data['Pos'] < 0)
    right_half_mask = (data['Pos'] >= 0)
    
    # Get run columns
    run_cols = slice(3, -4)
    runs = data.iloc[:, run_cols]
    
    # Initialize dictionary to store scaling factors
    scaling_factors = {}
    
    # Calculate scaling factors using flank regions for each strand
    for strand in data['Strand'].unique():
        strand_mask = (data['Strand'] == strand)
        
        # Left flank scaling (calculated from flanks)
        left_data = data[strand_mask & left_flank_mask]
        if not left_data.empty:
            left_factor = (left_data['Dmg'].mean() / 
                         left_data['mean'].mean())
            scaling_factors[f'{strand}_left'] = left_factor
            
        # Right flank scaling (calculated from flanks)
        right_data = data[strand_mask & right_flank_mask]
        if not right_data.empty:
            right_factor = (right_data['Dmg'].mean() / 
                          right_data['mean'].mean())
            scaling_factors[f'{strand}_right'] = right_factor
    
    print("Scaling Factors:", scaling_factors)
    
    # Apply scaling factors to entire left/right halves
    for strand in data['Strand'].unique():
        strand_mask = (data['Strand'] == strand)
        
        # Scale entire left half
        left_scaling_mask = strand_mask & left_half_mask
        if left_scaling_mask.any():
            data.loc[left_scaling_mask, data.columns[run_cols]] *= scaling_factors[f'{strand}_left']
        
        # Scale entire right half
        right_scaling_mask = strand_mask & right_half_mask
        if right_scaling_mask.any():
            data.loc[right_scaling_mask, data.columns[run_cols]] *= scaling_factors[f'{strand}_right']
    
    # Get updated runs view
    runs = data.iloc[:, run_cols]
    exp = data['Dmg']
    
    # Initialize p-value arrays
    p_values_top = np.zeros(len(data))
    p_values_bottom = np.zeros(len(data))
    
    # Calculate p-values position by position
    for i in range(len(data)):
        observed = exp.iloc[i]
        simulated = runs.iloc[i].values
        
        # If all values are identical (including observed), set p-values to 1
        if np.allclose(simulated, observed, rtol=1e-5, atol=1e-5):
            p_values_top[i] = 1.0
            p_values_bottom[i] = 1.0
        else:
            # Calculate empirical p-values only when there's actual variation
            p_values_top[i] = (simulated > observed).mean()
            p_values_bottom[i] = (simulated < observed).mean()
    
    # Apply FDR correction
    data['P-Value-Top'] = smstats.fdrcorrection(p_values_top)[1]
    data['P-Value-Bottom'] = smstats.fdrcorrection(p_values_bottom)[1]
    
    # Recalculate all summary statistics after scaling
    data['min'] = runs.min(axis=1)
    data['max'] = runs.max(axis=1)
    data['median'] = runs.median(axis=1)
    data['mean'] = runs.mean(axis=1)
    data['std'] = runs.std(axis=1)
    data['z'] = (data['Dmg'] - data['mean']) / (data['std'] + 1e-9)


def add_stats_muts(data, num_runs):
    data.reset_index(inplace=True) 
    
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)

    # Create mask for flanks
    flank_mask = ((data['Pos'] < motif_range[0] - 5) | (data['Pos'] > motif_range[1] + 5))
    
    # Get run columns with correct slice
    run_cols = slice(2, -4)  # Changed from -3 to -4 to account for all summary statistics
    
    if flank_mask.any():
        # Calculate factor using masked data
        factor = (data.loc[flank_mask, 'Mut'].mean() / data.loc[flank_mask, 'mean'].mean())
        print(f'Scaling Factor: {factor}')
        
        # Modify run columns in place
        data.iloc[:, run_cols] *= factor
    
    # Get runs view after modification
    runs = data.iloc[:, run_cols]
    exp = data['Mut']
        
    # Calculate p-values
    bool_comparison_top = runs.gt(exp, axis=0)
    bool_comparison_bottom = runs.lt(exp, axis=0)

    # Calculate FDR-corrected p-values
    data['P-Value-Top'] = smstats.fdrcorrection(
        bool_comparison_top.sum(axis=1) / num_runs)[1]
    data['P-Value-Bottom'] = smstats.fdrcorrection(
        bool_comparison_bottom.sum(axis=1) / num_runs)[1]

    # Recalculate all summary statistics after scaling
    data['min'] = runs.min(axis=1)
    data['max'] = runs.max(axis=1)
    data['median'] = runs.median(axis=1)
    data['mean'] = runs.mean(axis=1)
    data['std'] = runs.std(axis=1)
    data['z'] = (data['Mut'] - data['mean']) / (data['std'] + 1e-9)

def add_stats_mut_split(data, num_runs):
    data.reset_index(inplace=True)
    
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)
    
    # Create masks for flanks (for calculating scaling factors)
    left_flank_mask = (data['Pos'] < motif_range[0] - 5)
    right_flank_mask = (data['Pos'] > motif_range[1] + 5)
    
    # Create masks for applying scaling (entire left/right halves)
    left_half_mask = (data['Pos'] < 0)
    right_half_mask = (data['Pos'] >= 0)
    
    # Get run columns
    run_cols = slice(2, -4)
    runs = data.iloc[:, run_cols]
    
    # Initialize dictionary to store scaling factors
    scaling_factors = {}
    
    # Calculate left scaling factor from left flanks
    left_data = data[left_flank_mask]
    if not left_data.empty:
        left_factor = (left_data['Mut'].mean() / 
                      left_data['mean'].mean())
        scaling_factors['left'] = left_factor
    
    # Calculate right scaling factor from right flanks
    right_data = data[right_flank_mask]
    if not right_data.empty:
        right_factor = (right_data['Mut'].mean() / 
                       right_data['mean'].mean())
        scaling_factors['right'] = right_factor
    
    print("Scaling Factors:", scaling_factors)
    
    # Apply scaling factors to entire halves
    # Scale left half
    left_scaling_mask = left_half_mask
    if left_scaling_mask.any():
        data.loc[left_scaling_mask, data.columns[run_cols]] *= scaling_factors['left']
    
    # Scale right half
    right_scaling_mask = right_half_mask
    if right_scaling_mask.any():
        data.loc[right_scaling_mask, data.columns[run_cols]] *= scaling_factors['right']
    
    # Get updated runs view
    runs = data.iloc[:, run_cols]
    exp = data['Mut']
    
    # Initialize p-value arrays
    p_values_top = np.zeros(len(data))
    p_values_bottom = np.zeros(len(data))
    
    # Calculate p-values position by position
    for i in range(len(data)):
        observed = exp.iloc[i]
        simulated = runs.iloc[i].values
        
        # If all values are identical (including observed), set p-values to 1
        if np.allclose(simulated, observed, rtol=1e-5, atol=1e-5):
            p_values_top[i] = 1.0
            p_values_bottom[i] = 1.0
        else:
            # Calculate empirical p-values only when there's actual variation
            p_values_top[i] = (simulated > observed).mean()
            p_values_bottom[i] = (simulated < observed).mean()
    
    # Apply FDR correction
    data['P-Value-Top'] = smstats.fdrcorrection(p_values_top)[1]
    data['P-Value-Bottom'] = smstats.fdrcorrection(p_values_bottom)[1]
    
    # Recalculate all summary statistics after scaling
    data['min'] = runs.min(axis=1)
    data['max'] = runs.max(axis=1)
    data['median'] = runs.median(axis=1)
    data['mean'] = runs.mean(axis=1)
    data['std'] = runs.std(axis=1)
    data['z'] = (data['Mut'] - data['mean']) / (data['std'] + 1e-9)
    
    return scaling_factors

def scale_mut_mutagenic(data, proj_col, act_col):

    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)

    flanks = data[((data.index < motif_range[0] - 5) | (data.index > motif_range[1] + 5))]
    if len(flanks) > 0:
        factor = flanks[act_col].mean() / flanks[proj_col].mean()
        data[proj_col] = data[proj_col] * factor

        #data['Mut'] = data['Mut'] * flanks.mean()
        print(f'Scaling Factor: {factor}')
    
def plot_data_with_stats(data, data_cols, tf_name, window, ax, title_prefix="", y_label="Counts"):
    """Generic plotting function for damage/mutation data with statistics
    
    Args:
        data: DataFrame containing the data to plot
        data_cols: List of column names to plot (e.g. ['Dmg', 'Mut'])
        tf_name: Name of transcription factor
        window: Window size
        ax: Matplotlib axis to plot on
        title_prefix: Prefix for plot title
        y_label: Label for y-axis
    """
    high = 0
    
    # Plot each data column
    for col in data_cols:
        if 'Strand' in data.columns:
            for strand, group in data.groupby('Strand'):
                label = "Opposite" if strand == 'Diff' else 'Same'
                ax.plot(group['Pos'], group[col], label=label)
                ax.scatter(group['Pos'], group[col], label=None)
                
                if 'median' in group:
                    ax.plot(group['Pos'], group['median'], linestyle='dashed', alpha=0.5)
                if 'min' in group and 'max' in group:
                    ax.fill_between(group['Pos'], group['min'], group['max'], alpha=0.5, linestyle='dashed')
                    high = max(high, group['max'].max())
                high = max(high, group[col].max())
        else:
            ax.plot(data['Pos'], data[col], label=col)
            ax.scatter(data['Pos'], data[col], label=None)
            if 'median' in data:
                ax.plot(data['Pos'], data['median'], linestyle='dashed', alpha=0.5)
            if 'min' in data and 'max' in data:
                ax.fill_between(data['Pos'], data['min'], data['max'], alpha=0.5, linestyle='dashed')
                high = max(high, data['max'].max())
            high = max(high, data[col].max())

    # Set plot formatting
    ax.set_xlim(-15.5, 15.5)
    max_y = (high // 500 + 1) * 500
    ax.set_ylim(-250, max_y)
    
    ax.set_xlabel('Position')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.set_ylabel(y_label)
    ax.set_title(f'{title_prefix} Profile for {tf_name}')

    # Add motif region shading and grid lines
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)
    ax.axvspan(motif_range[0], motif_range[1], alpha=0.35, color='grey')

    for x in np.arange(-15, 15, 1):
        ax.axvline(x, linestyle='dashed', color='grey', linewidth=0.5)
    for y in np.arange(-250, max_y, 250):
        ax.axhline(y, linestyle='dashed', color='grey', linewidth=0.5)

    # Add motif boundary lines
    ax.axvline(motif_range[0], linestyle='solid', color='#2b2b2b', linewidth=1)
    ax.axvline(motif_range[1], linestyle='solid', color='#2b2b2b', linewidth=1)
    ax.axvline(0, linestyle='dashed', color='#2b2b2b', linewidth=1)

    ax.legend(loc='upper right')

def generate_logo(res_cur, plus_counts, minus_counts, tf_name, ax):
    """Generate sequence logo plot with validation"""
    # Debug information
    print(f"Number of plus sequences: {len(plus_counts)}")
    print(f"Number of minus sequences: {len(minus_counts)}")
    
    # Validate sequence lengths
    plus_lengths = [len(seq) for seq in plus_counts]
    minus_lengths = [len(seq) for seq in minus_counts]
    
    if len(set(plus_lengths)) > 1:
        print("Warning: Inconsistent plus sequence lengths:", set(plus_lengths))
        problematic = [(i, len(s), s) for i, s in enumerate(plus_counts) if len(s) != 31]
        print("Problematic plus sequences:", problematic[:5])  # Show first 5 issues
        
    if len(set(minus_lengths)) > 1:
        print("Warning: Inconsistent minus sequence lengths:", set(minus_lengths))
        problematic = [(i, len(s), s) for i, s in enumerate(minus_counts) if len(s) != 31]
        print("Problematic minus sequences:", problematic[:5])  # Show first 5 issues
    
    # Filter out invalid sequences
    plus_counts = [seq for seq in plus_counts if len(seq) == 31]
    minus_counts = [seq for seq in minus_counts if len(seq) == 31]
    
    if not plus_counts or not minus_counts:
        raise ValueError("No valid sequences after filtering")
    
    pos_plus = [i for i in range(-15, 16)]
    plus_motif = lm.alignment_to_matrix(sequences=plus_counts, to_type='counts', characters_to_ignore='.-X')
    plus_motif['pos'] = pos_plus
    plus_motif.set_index('pos', inplace=True)

    pos_minus = [i for i in range(15, -16, -1)]
    minus_motif = lm.alignment_to_matrix(sequences=minus_counts, to_type='counts', characters_to_ignore='.-X')
    minus_motif['pos'] = pos_minus
    minus_motif.set_index('pos', inplace=True)

    comb = pd.DataFrame()
    comb['pos'] = pos_plus
    comb.set_index('pos', inplace=True)

    for nuc in ['A', 'C', 'G', 'T']:
        comb[nuc] = plus_motif[nuc] - minus_motif[nuc]

    logo = lm.Logo(df=comb, fade_below=0, shade_below=0, ax=ax)
    
    logo.ax.set_yticklabels([])
    logo.ax.axes.get_xaxis().set_visible(False)
    logo.ax.set_ylabel("Empirical Binding Motif", labelpad=12)
    logo.ax.axvspan(motif_range[0] - 0.5, motif_range[1] + 0.5, alpha=0.35, color='grey')

    for x in np.arange(-15, 15, 1):
        logo.ax.axvline(x, linestyle='dashed', color='grey', linewidth=0.5)

    ax.axvline(motif_range[0] - 0.5, linestyle='solid', color='#2b2b2b', linewidth=1)
    ax.axvline(motif_range[1] + 0.5, linestyle='solid', color='#2b2b2b', linewidth=1)
    ax.axvline(0, linestyle='dashed', color='#2b2b2b', linewidth=1)

def plot_data_with_stats_adaptive(data, data_cols, tf_name, window, ax, title_prefix="", y_label="Counts"):
    """Generic plotting function for damage/mutation data with statistics"""
    high = 0
    
    # Define colors 
    strand_colors = {'Same': '#1f77b4', 'Diff': '#ff7f0e'}  # For damage case
    single_color = '#1f77b4'  # For mutation case - using the same blue as strand_colors['Same']
    
    # Plot each data column
    for col in data_cols:
        if 'Strand' in data.columns:
            # Stranded case (damage)
            for strand, group in data.groupby('Strand'):
                label = "Opposite Strand" if strand == 'Diff' else 'Motif Strand'
                color = strand_colors[strand]
                ax.plot(group['Pos'], group[col], label=label, color=color)
                
                # Determine marker for each point based on p-values
                for idx, row in group.iterrows():
                    marker = 'o'
                    markersize = 8
                    if row['P-Value-Top'] < 0.05:
                        marker = '^'
                        markersize = 10
                    elif row['P-Value-Bottom'] < 0.05:
                        marker = 'v'
                        markersize = 10
                    ax.scatter(row['Pos'], row[col], marker=marker, s=markersize**2, 
                             color=color, label=None)
                
                if 'median' in group:
                    ax.plot(group['Pos'], group['median'], linestyle='dashed', 
                           alpha=0.5, color=color)
                if 'min' in group and 'max' in group:
                    ax.fill_between(group['Pos'], group['min'], group['max'], 
                                  alpha=0.5, color=color, linestyle='dashed')
                    high = max(high, group['max'].max())
                high = max(high, group[col].max())
        else:
            # Unstranded case (mutation)
            ax.plot(data['Pos'], data[col], label=col, color=single_color)
            
            # Determine marker for each point based on p-values
            for idx, row in data.iterrows():
                marker = 'o'
                markersize = 8
                if row['P-Value-Top'] < 0.05:
                    marker = '^'
                    markersize = 10
                elif row['P-Value-Bottom'] < 0.05:
                    marker = 'v'
                    markersize = 10
                ax.scatter(row['Pos'], row[col], marker=marker, s=markersize**2, 
                         color=single_color, label=None)
            
            if 'median' in data:
                ax.plot(data['Pos'], data['median'], linestyle='dashed', 
                       alpha=0.5, color=single_color)
            if 'min' in data and 'max' in data:
                ax.fill_between(data['Pos'], data['min'], data['max'], 
                              alpha=0.5, color=single_color, linestyle='dashed')
                high = max(high, data['max'].max())
            high = max(high, data[col].max())

    # Rest of the formatting code remains the same
    ax.set_xlim(-15.5, 15.5)
    max_y = (high // 500 + 1) * 500
    ax.set_ylim(-250, max_y)
    
    ax.set_xlabel('Position')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.set_ylabel(y_label)
    ax.set_title(f'{title_prefix} Profile for {tf_name}')

    # Add motif region shading and grid lines
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)
    ax.axvspan(motif_range[0], motif_range[1], alpha=0.35, color='grey')

    for x in np.arange(-15, 15, 1):
        ax.axvline(x, linestyle='dashed', color='grey', linewidth=0.5)
    for y in np.arange(-250, max_y, 250):
        ax.axhline(y, linestyle='dashed', color='grey', linewidth=0.5)

    # Add motif boundary lines
    ax.axvline(motif_range[0], linestyle='solid', color='#2b2b2b', linewidth=1)
    ax.axvline(motif_range[1], linestyle='solid', color='#2b2b2b', linewidth=1)
    ax.axvline(0, linestyle='dashed', color='#2b2b2b', linewidth=1)

    ax.legend(loc='upper right')

def plot_pvalue_zscore_logo(data, tf_name, plus_counts, minus_counts, window):
    """Generate 3-panel plot with p-values, z-scores, and sequence logo"""
    plt.rcParams['font.size'] = 36
    fig, axs = plt.subplots(nrows=3, ncols=1, layout='constrained', figsize=(18, 24))
    
    # Define colors
    strand_colors = {'Same': '#1f77b4', 'Diff': '#ff7f0e'}
    single_color = '#1f77b4'
    
    # Plot 1: -log10 p-values
    if 'Strand' in data.columns:
        # Stranded case (damage)
        for strand, group in data.groupby('Strand'):
            label = "Opposite Strand" if strand == 'Diff' else 'Motif Strand'
            color = strand_colors[strand]
            
            neg_log_p_top = -np.log10(group['P-Value-Top'] + 0.00001)
            neg_log_p_bottom = -np.log10(group['P-Value-Bottom'] + 0.00001)
            
            axs[0].scatter(group['Pos'], neg_log_p_top, marker='^', s=100, 
                          color=color, label=f"{label} (Enriched)")
            axs[0].scatter(group['Pos'], neg_log_p_bottom, marker='v', s=100, 
                          color=color, label=f"{label} (Depleted)")
    else:
        # Unstranded case (mutation)
        neg_log_p_top = -np.log10(data['P-Value-Top'] + 0.00001)
        neg_log_p_bottom = -np.log10(data['P-Value-Bottom'] + 0.00001)
        
        axs[0].scatter(data['Pos'], neg_log_p_top, marker='^', s=100, 
                      color=single_color, label="Enriched")
        axs[0].scatter(data['Pos'], neg_log_p_bottom, marker='v', s=100, 
                      color=single_color, label="Depleted")
    
    # Add significance threshold line
    sig_threshold = -np.log10(0.05)
    axs[0].axhline(y=sig_threshold, color='red', linestyle='--', label='p = 0.05')
    axs[0].axhline(y=0, color='#2b2b2b', linestyle='-', linewidth=1)
    axs[0].set_ylabel("-log10(p-value)")
    axs[0].legend(loc='upper right')
    
    # Plot 2: Z-scores
    if 'Strand' in data.columns:
        # Stranded case (damage)
        for strand, group in data.groupby('Strand'):
            label = "Opposite Strand" if strand == 'Diff' else 'Motif Strand'
            color = strand_colors[strand]
            axs[1].scatter(group['Pos'], group['z'], s=100, color=color, label=label)
            axs[1].plot(group['Pos'], group['z'], color=color, alpha=0.5)
    else:
        # Unstranded case (mutation)
        axs[1].scatter(data['Pos'], data['z'], s=100, color=single_color, label='Z-score')
        axs[1].plot(data['Pos'], data['z'], color=single_color, alpha=0.5)
    
    axs[1].axhline(y=0, color='#2b2b2b', linestyle='-', linewidth=1)
    axs[1].set_ylabel("Z-score")
    axs[1].legend(loc='upper right')
    
    # Plot 3: Sequence logo
    generate_logo(data, plus_counts, minus_counts, tf_name, axs[2])
    
    # Add shared formatting
    for ax in axs:
        ax.set_xlim(-15.5, 15.5)
        ax.grid(True, linestyle='--', alpha=0.3)
        
        # Add motif region shading and boundary lines
        global motif_range
        motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)
        ax.axvspan(motif_range[0], motif_range[1], alpha=0.35, color='grey')
        ax.axvline(motif_range[0], linestyle='solid', color='#2b2b2b', linewidth=1)
        ax.axvline(motif_range[1], linestyle='solid', color='#2b2b2b', linewidth=1)
        ax.axvline(0, linestyle='dashed', color='#2b2b2b', linewidth=1)
        
        for x in np.arange(-15, 15, 1):
            ax.axvline(x, linestyle='dashed', color='grey', linewidth=0.5)
    
    axs[2].set_xlabel('Position')
    plt.show()
    fig.savefig(f'../../output/uv_rescaled_mut/pvalue_zscore_logo_{tf_name}_{window}.png', dpi=600)

def plot_stacked_analysis(data, data_cols, tf_name, window, plus_counts, minus_counts, mode='damage', title_prefix="", save_path="../../output/"):
    """Generate a stacked figure with damage/mutation data, p-values, z-scores, and sequence logo"""
    plt.rcParams['font.size'] = 42  # Base font size
    
    # Process TF name to replace _number with cnumber
    if '_' in tf_name and tf_name.split('_')[-1].isdigit():
        base_name = '_'.join(tf_name.split('_')[:-1])
        number = tf_name.split('_')[-1]
        tf_name = f"{base_name}c{number}"
    
    # Create figure with updated dimensions and height ratios
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(30, 21), 
                           gridspec_kw={'height_ratios': [7, 6, 5, 3]},
                           layout='constrained')
    
    # Define colors
    strand_colors = {'Same': '#1f77b4', 'Diff': '#ff7f0e'}
    single_color = '#1f77b4'
    
    # Plot 1: Damage/Mutation Data
    high = 0
    for col in data_cols:
        if 'Strand' in data.columns:
            for strand, group in data.groupby('Strand'):
                label = "Opposite Strand" if strand == 'Diff' else 'Motif Strand'
                color = strand_colors[strand]
                axs[0].plot(group['Pos'], group[col], label=label, color=color)
                
                # Larger markers for significance
                for idx, row in group.iterrows():
                    marker = 'o'
                    markersize = 20 # Increased base size
                    if row['P-Value-Top'] < 0.05:
                        marker = '^'
                        markersize = 20  # Increased arrow size
                    elif row['P-Value-Bottom'] < 0.05:
                        marker = 'v'
                        markersize = 20
                    axs[0].scatter(row['Pos'], row[col], marker=marker, s=markersize**2, 
                                 color=color, label=None)
                
                if 'median' in group:
                    axs[0].plot(group['Pos'], group['median'], linestyle='dashed', 
                              alpha=0.5, color=color)
                if 'min' in group and 'max' in group:
                    axs[0].fill_between(group['Pos'], group['min'], group['max'], 
                                      alpha=0.5, color=color, linestyle='dashed')
                    high = max(high, group['max'].max())
                high = max(high, group[col].max())
        else:
            # Unstranded case (similar structure with single color)
            axs[0].plot(data['Pos'], data[col], label=col, color=single_color)
            for idx, row in data.iterrows():
                marker = 'o'
                markersize = 20
                if row['P-Value-Top'] < 0.05:
                    marker = '^'
                    markersize = 20
                elif row['P-Value-Bottom'] < 0.05:
                    marker = 'v'
                    markersize = 20
                axs[0].scatter(row['Pos'], row[col], marker=marker, s=markersize**2, 
                             color=single_color, label=None)
            
            if 'median' in data:
                axs[0].plot(data['Pos'], data['median'], linestyle='dashed', 
                           alpha=0.5, color=single_color)
            if 'min' in data and 'max' in data:
                axs[0].fill_between(data['Pos'], data['min'], data['max'], 
                                  alpha=0.5, color=single_color, linestyle='dashed')
                high = max(high, data['max'].max())
            high = max(high, data[col].max())

    # Plot 2: -log10 p-values
    if 'Strand' in data.columns:
        for strand, group in data.groupby('Strand'):
            label = "Opposite Strand" if strand == 'Diff' else 'Motif Strand'
            color = strand_colors[strand]
            
            neg_log_p_top = -np.log10(group['P-Value-Top'] + 0.00001)
            neg_log_p_bottom = -np.log10(group['P-Value-Bottom'] + 0.00001)
            
            axs[1].scatter(group['Pos'], neg_log_p_top, marker='^', s=225,  # Increased marker size
                          color=color, label=f"{label} (Enriched)")
            axs[1].scatter(group['Pos'], neg_log_p_bottom, marker='v', s=225,
                          color=color, label=f"{label} (Depleted)")
    else:
        neg_log_p_top = -np.log10(data['P-Value-Top'] + 0.00001)
        neg_log_p_bottom = -np.log10(data['P-Value-Bottom'] + 0.00001)
        
        axs[1].scatter(data['Pos'], neg_log_p_top, marker='^', s=225,
                      color=single_color, label="Enriched")
        axs[1].scatter(data['Pos'], neg_log_p_bottom, marker='v', s=225,
                      color=single_color, label="Depleted")

    # Plot 3: Z-scores
    if 'Strand' in data.columns:
        for strand, group in data.groupby('Strand'):
            label = "Opposite Strand" if strand == 'Diff' else 'Motif Strand'
            color = strand_colors[strand]
            axs[2].scatter(group['Pos'], group['z'], s=225, color=color, label=label)
            axs[2].plot(group['Pos'], group['z'], color=color, alpha=0.5)
    else:
        axs[2].scatter(data['Pos'], data['z'], s=225, color=single_color, label='Z-score')
        axs[2].plot(data['Pos'], data['z'], color=single_color, alpha=0.5)

    # Plot 4: Sequence logo
    generate_logo(data, plus_counts, minus_counts, tf_name, axs[3])

    # Set titles and labels
    max_y = (high // 500 + 1) * 500
    axs[0].set_ylim(-250, max_y)
    axs[0].set_ylabel("Count")  # Simplified label
    axs[0].set_title(f"{title_prefix} Profile for {tf_name}")
    
    # Add significance threshold line to p-value plot
    sig_threshold = -np.log10(0.05)
    axs[1].axhline(y=sig_threshold, color='red', linestyle='--', label='q = 0.05')
    axs[1].set_ylabel("-log10(q-value)")
    
    axs[2].set_ylabel("Z-score")
    axs[2].axhline(y=0, color='#2b2b2b', linestyle='-', linewidth=1)

    axs[3].set_xlabel("Position")
    
    # Remove the empirical binding motif label
    axs[3].set_ylabel("")

    # Shared formatting for all subplots with larger legend
    global motif_range
    motif_range = (-(math.ceil(motif_length / 2) - 1), motif_length // 2)

    for i, ax in enumerate(axs):
        ax.set_xlim(-15.5, 15.5)
        ax.grid(True, linestyle='--', alpha=0.3)
        
        if i == 3:  # Logo plot (last subplot)
            # Add motif region shading and boundary lines with offset for logo only
            ax.axvspan(motif_range[0] - 0.5, motif_range[1] + 0.5, alpha=0.35, color='grey')
            ax.axvline(motif_range[0] - 0.5, linestyle='solid', color='#2b2b2b', linewidth=1)
            ax.axvline(motif_range[1] + 0.5, linestyle='solid', color='#2b2b2b', linewidth=1)
        else:
            # Regular motif region shading and boundary lines for other plots
            ax.axvspan(motif_range[0], motif_range[1], alpha=0.35, color='grey')
            ax.axvline(motif_range[0], linestyle='solid', color='#2b2b2b', linewidth=1)
            ax.axvline(motif_range[1], linestyle='solid', color='#2b2b2b', linewidth=1)
        
        ax.axvline(0, linestyle='dashed', color='#2b2b2b', linewidth=1)
        
        for x in np.arange(-15, 15, 1):
            ax.axvline(x, linestyle='dashed', color='grey', linewidth=0.5)
        
        ax.legend(loc='upper right', fontsize=32)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{save_path}/{mode}_complete_analysis_{tf_name}_{window}.png', dpi=600, bbox_inches='tight')

def summary_visuals(results_damage, results_mutation, tf_name, window, mode='damage'):
    """Generate the complete stacked visualization"""
    if results_damage is not None:
        plot_stacked_analysis(results_damage, ['Dmg'], tf_name, window, 
                            plus_tf, minus_tf, mode='damage',
                            title_prefix="Damage", 
                            save_path="../../output/uv_rescaled_dmg")
    else:
        plot_stacked_analysis(results_mutation, ['Mut'], tf_name, window,
                            plus_tf, minus_tf, mode='mutation',
                            title_prefix="Mutation",
                            save_path="../../output/uv_rescaled_mut")

def write_results(tf_name, results, alpha=0.05):
    from pathlib import Path

    # Filter data once and reuse
    relevant = results[(results["Pos"] > -15) & (results["Pos"] < 15)][["Pos", "P-Value-Top", "P-Value-Bottom", "Strand", "z"]]
    
    # Calculate stats
    min_pval_top = relevant['P-Value-Top'].min()
    significant_top = relevant[relevant['P-Value-Top'] < alpha]
    num_sig_top = significant_top.shape[0]

    min_pval_bottom = relevant['P-Value-Bottom'].min() 
    significant_bottom = relevant[relevant['P-Value-Bottom'] < alpha]
    num_sig_bottom = significant_bottom.shape[0]

    max_z = relevant['z'].max()
    min_z = relevant['z'].min()

    # Write results using context manager
    output_path = Path('../../output/uv_rescaled_dmg') / f'{tf_name}.out'
    with open(output_path, 'w') as out:
        out.write(f'SIGNIFICANCE ANALYSIS FOR {tf_name}\n\n')
        out.write('Min-Top\t\tNum-Top\t\tMin-Bottom\t\tNum-Bottom\t\tLargest Z\t\tSmallest Z\n')
        out.write(f'{min_pval_top}\t\t{num_sig_top}\t\t{min_pval_bottom}\t\t{num_sig_bottom}\t\t{max_z}\t\t{min_z}\n')
        
        out.write('\nSIGNIFICANT POSITIONS - ENRICHMENT\n')
        for _, row in significant_top.iterrows():
            out.write(f"{row['Pos']}\t{row['P-Value-Top']}\t{row['Strand']}\t{row['z']}\n")

        out.write('\nSIGNIFICANT POSITIONS - PROTECTION\n') 
        for _, row in significant_bottom.iterrows():
            out.write(f"{row['Pos']}\t{row['P-Value-Bottom']}\t{row['Strand']}\t{row['z']}\n")

def write_results_muts(tf_name, results, alpha=0.05):
    from pathlib import Path

    # Filter data once and reuse - removed Strand from selection
    relevant = results[(results["Pos"] > -15) & (results["Pos"] < 15)][["Pos", "P-Value-Top", "P-Value-Bottom", "z"]]
    
    # Calculate stats
    min_pval_top = relevant['P-Value-Top'].min()
    significant_top = relevant[relevant['P-Value-Top'] < alpha]
    num_sig_top = significant_top.shape[0]

    min_pval_bottom = relevant['P-Value-Bottom'].min() 
    significant_bottom = relevant[relevant['P-Value-Bottom'] < alpha]
    num_sig_bottom = significant_bottom.shape[0]

    max_z = relevant['z'].max()
    min_z = relevant['z'].min()

    # Write results using context manager
    output_path = Path('../../output/uv_rescaled_mut/') / f'{tf_name}.out'
    with open(output_path, 'w') as out:
        out.write(f'SIGNIFICANCE ANALYSIS FOR {tf_name} MUTATIONS\n\n')
        out.write('Min-Top\t\tNum-Top\t\tMin-Bottom\t\tNum-Bottom\t\tLargest Z\t\tSmallest Z\n')
        out.write(f'{min_pval_top}\t\t{num_sig_top}\t\t{min_pval_bottom}\t\t{num_sig_bottom}\t\t{max_z}\t\t{min_z}\n')
        
        out.write('\nSIGNIFICANT POSITIONS - INCREASED MUTATIONS\n')
        for _, row in significant_top.iterrows():
            out.write(f"{row['Pos']}\t{row['P-Value-Top']}\t{row['z']}\n")  # Removed Strand

        out.write('\nSIGNIFICANT POSITIONS - DECREASED MUTATIONS\n')
        for _, row in significant_bottom.iterrows():
            out.write(f"{row['Pos']}\t{row['P-Value-Bottom']}\t{row['z']}\n")  # Removed Strand

def compare_mutations(counts, mutagenic, actual, tf_name, alpha=0.05, alpha_2=0.5):

    p_combined = pd.concat([mutagenic["Pos"], mutagenic["Mut"], mutagenic["z"], mutagenic["P-Value-Top"], mutagenic["P-Value-Bottom"], actual["Mut"], actual["z"], actual["P-Value-Top"], actual["P-Value-Bottom"]], axis=1)
    p_combined.columns = ['Pos', 'Counts-Potential', 'Z-Potential', 'PVT-Potential', 'PVB-Potential', 'Counts-Actual', 'Z-Actual', 'PVT-Actual', 'PVB-Actual']
    scale_mut_mutagenic(p_combined, 'Counts-Potential', 'Counts-Actual')
    
    significant = p_combined[(p_combined["Pos"] > -15) & (p_combined["Pos"] < 15)]
    
    significant.to_csv(f'stats/atac_mut/{tf_name}_summary.csv')

    pvt_actual, pvt_potential = significant['PVT-Actual'] < alpha, significant['PVT-Potential'] < alpha
    pvb_actual, pvb_potential = significant['PVB-Actual'] < alpha, significant['PVB-Potential'] < alpha

    pvt_insig_actual, pvt_insig_potential = significant['PVT-Actual'] > alpha_2, significant['PVT-Potential'] > alpha_2
    pvb_insig_actual, pvb_insig_potential = significant['PVB-Actual'] > alpha_2, significant['PVB-Potential'] > alpha_2

    # Print tab-delimited table
    print("Test\tBoth Sig.\tOnly Potential\tOnly Actual\tNeither")
    print(f"Enrichment\t{sum(pvt_actual & pvt_potential)}\t{sum(pvt_insig_actual & pvt_potential)}\t{sum(pvt_actual & pvt_insig_potential)}\t{sum(pvt_insig_actual & pvt_insig_potential)}")
    print(f"Depletion\t{sum(pvb_actual & pvb_potential)}\t{sum(pvb_insig_actual & pvb_potential)}\t{sum(pvb_actual & pvb_insig_potential)}\t{sum(pvb_insig_actual & pvb_insig_potential)}")
    print("")
    # print("\n")
    # print(significant)

def permutation_statistic(x, y):
    dof = len(x) - 2
    rs = stats.spearmanr(x, y).statistic
    transformed = rs * np.sqrt(dof / ((rs+1.0)*(1.0-rs)))
    return transformed

def analyze(tf_name, num_runs, run_id, window=20):
    global plus_tf
    plus_tf = []

    global minus_tf
    minus_tf = []

    #motif_length = 0

    #tree = process_intervals(read_tf_file("tf_exp/accessible_ETS1_plus.bed") + read_tf_file("tf_exp/accessible_ETS1_minus.bed"), 300)
    # tree = process_intervals(read_tf_file(f'/home/users/bc301/scripts/intervals/accessibility/processed_clusters/a549/{tf_name}.bed'), window)
    tree = process_intervals(read_tf_file(f'../../data/tfbs/fib_atac/{tf_name}/{tf_name}_plus_high.bed') + read_tf_file(f'../../data/tfbs/fib_atac/{tf_name}/{tf_name}_minus_high.bed'), window)
    # tree = process_intervals(read_tf_file(f'../../data/tfbs/a549/{tf_name}/{tf_name}_plus_high.bed') + read_tf_file(f'../../data/tfbs/a549/{tf_name}/{tf_name}_minus_high.bed'), window)
    # tree = process_intervals(read_tf_file(f'ETS_1_pyr3.bed'), window) # test wyrick fig. 4
    # intervals = hold_interval_indices('acc_regions/atac_accessibility_sorted.bed')
    # intervals = hold_interval_indices('/usr/xtmp/bc301/hana_data/for_bo/WT_CSB_hg19_idr_conservative_summits_150bp.bed')

    # FOR DAMAGE

    # intervals = hold_interval_indices('/home/users/bc301/scripts/intervals/accessibility/accessibility_sorted.bed')
    
    intervals = hold_interval_indices('../../data/raw/atac_150bp_intervals_merged.bed')
    # intervals = hold_interval_indices('../../data/raw/a549_regions_merged.bed')

    # DAMAGE ANALYSIS !!!
    total_runs = build_runs(num_runs, intervals, tree, window)
    intersection_dmgs = process_dmgs('../../data/damages/atac_nb_plus.bed', '../../data/damages/atac_nb_minus.bed', tree, window) #= 
    # intersection_dmgs = process_dmgs('../../data/damages/a549_bpde_1_2_plus.bed', '../../data/damages/a549_bpde_1_2_minus.bed', tree, window) #= 

    combined_dmg = pd.concat([intersection_dmgs, total_runs], axis=1)

    # #add_stats(combined_dmg, num_runs)
    add_stats_split(combined_dmg, num_runs)
    combined_dmg.to_csv(f"../../output/uv_rescaled_dmg/{tf_name}_{num_runs}_{window}_dmg_raw_{run_id}_nb.csv")
    write_results(tf_name, combined_dmg)
    summary_visuals(combined_dmg, None, tf_name, window, mode='damage')




    # intersection_dmgs = process_dmgs('damages/atac_accessible_naked_plus.bed', 'damages/atac_accessible_naked_minus.bed', tree, window) #= 
    # total_runs = build_runs(num_runs, intervals, tree, window)
    # combined_dmg = pd.concat([intersection_dmgs, total_runs], axis=1)
    # # combined_dmg.reset_index(inplace=True)
    # combined_dmg.rename(columns={'index': 'Pos'}, inplace=True)
    # add_stats(combined_dmg, num_runs)
    # combined_dmg.to_csv(f"../../output/preprocessed/raw/{tf_name}_{num_runs}_{window}_dmg_raw.csv")

    # FOR MUTAGENIC (POTENTIAL DAMAGE)
    # intersection_mutagenic = process_muts('mutations/mutagenic_combined_Conly.bed', tree, window)
    # total_mutagenic = build_muts(num_runs, intervals, tree, window, 'potential')

    # FOR ACTUAL MUTATIONS
    intersection_muts = process_muts('../../data/raw/atac_mutations_transitions_C_only.bed', tree, window)
    
    get_memory_usage("Before mutation build")
    total_muts = build_muts(num_runs, intervals, tree, window, 'actual')
    get_memory_usage("After mutation build")
    
    combined_mut = pd.concat([intersection_muts, total_muts], axis=1)
    add_stats_mut_split(combined_mut, num_runs)

    combined_mut.to_csv(f"../../output/uv_rescaled_mut/{tf_name}_{num_runs}_{window}_mutation.csv")
    summary_visuals(None, combined_mut, tf_name, window, mode='mutation')
    write_results_muts(tf_name, combined_mut)

    # comparison = pd.concat([intersection_mutagenic, intersection_muts], axis=1)
    # comparison.columns = ['Potential', 'Actual']
    # scale_mut_mutagenic(comparison, 'Potential', 'Actual')
    #comparison.to_csv(f"stats/mut/tf_testing/{tf_name}_{num_runs}_{window}_mutation_comparisons_a.csv")

    # combined = pd.concat([intersection_dmgs, total_runs], axis=1)
    # add_stats(combined, num_runs)

    #combined.to_csv(f"stats/mut/{tf_name}_{num_runs}_{window}_mutagenic_Conly.csv")

    # combined_dmg = pd.concat([intersection_mutagenic, total_mutagenic], axis=1)
    # add_stats_muts(combined_dmg, num_runs)


    
    # write_results(tf_name, combined_dmg)
    # summary_visuals(combined_dmg, None, tf_name, window, mode='damage')
    

    


    # # combined_dmg.to_csv(f"stats/mut/{tf_name}_{num_runs}_{window}_mutagenic_only_sim_a.csv")
    # # combined_mut.to_csv(f"stats/mut/{tf_name}_{num_runs}_{window}_mutation_Conly_sim_a.csv")
    # summary_visuals(None, combined_dmg, tf_name, window, mode='damage')
    # summary_visuals(None, combined_mut, tf_name, window, mode='mutation')

    # compare_mutations(comparison, combined_dmg, combined_mut, tf_name)

    # x_m, y_m = comparison['Potential'], comparison['Actual']
    # x_z, y_z = combined_dmg['z'], combined_mut['z']

    # rho_m, p_m = stats.spearmanr(x_m, y_m)
    # rho_z, p_z = stats.spearmanr(x_z, y_z)

    # print("\nWhole Statistics:")

    # print(f'Counts - Rho: {rho_m}, p-value: {p_m}')
    # print(f'Z-Score - Rho: {rho_z}, p-value: {p_z}')
    
    # ref_m = stats.permutation_test((x_m,y_m), permutation_statistic, alternative='greater', permutation_type='pairings')
    # ref_z = stats.permutation_test((x_z,y_z), permutation_statistic, alternative='greater', permutation_type='pairings')

    # print(f'Permutation Test - Counts: {ref_m.pvalue}, Z-Scores: {ref_z.pvalue}')

    

    # z_score_combined = pd.concat([combined_dmg['Pos'], combined_dmg['z'], combined_mut['z']], axis=1)
    # z_score_combined.columns = ['Pos', 'Potential Z', 'Actual Z']
    # #z_score_combined.to_csv(f"stats/mut/tf_testing/{tf_name}_{num_runs}_{window}_zscores_a.csv")

    # print("\nCore Statistics:")

    # comparison_core = comparison[(comparison.index >= motif_range[0]) & (comparison.index <= motif_range[1])]
    # z_score_core = z_score_combined[(z_score_combined['Pos'] >= motif_range[0]) & (z_score_combined['Pos'] <= motif_range[1])]
    
    # core_x_m, core_y_m = comparison_core['Potential'], comparison_core['Actual']
    # core_x_z, core_y_z = z_score_core['Potential Z'], z_score_core['Actual Z']

    # core_rho_m, core_p_m = stats.spearmanr(core_x_m, core_y_m)
    # core_rho_z, core_p_z = stats.spearmanr(core_x_z, core_y_z)

    # print(f'Counts - Rho: {core_rho_m}, p-value: {core_p_m}')
    # print(f'Z-Score - Rho: {core_rho_z}, p-value: {core_p_z}')

    # core_ref_m = stats.permutation_test((core_x_m,core_y_m), permutation_statistic, alternative='greater', permutation_type='pairings')
    # core_ref_z = stats.permutation_test((core_x_z,core_y_z), permutation_statistic, alternative='greater', permutation_type='pairings')

    # print(f'Permutation Test - Counts: {core_ref_m.pvalue}, Z-Scores: {core_ref_z.pvalue}')




    #summary_visuals(intersection_muts, combined, tf_name, window)
    
    
    
    # write_results(tf_name=tf_name, results=combined)

    

    
    #summary_visuals(combined, tf_name, window)

    
    

    

if __name__ == "__main__":
    # process_tfbs_file("ETS_1_plus.bed", "compressed_ETS1_plus.bed")
    # process_tfbs_file("ETS_1_minus.bed", "compressed_ETS1_minus.bed")
    #match_intervals(read_tf_file("tf_exp/accessible_ETS1_plus.bed"), get_intervals('accessibility_extended.bed'), "tf_exp/coded_ETS1_plus.bed")
    #match_intervals(read_dmg_file("accessible_plus.bed"), get_intervals('accessibility_sorted.bed'), "tf_exp/coded_dmg_plus.bed")
    #match_intervals(read_dmg_file("accessible_minus.bed"), get_intervals('accessibility_sorted.bed'), "tf_exp/coded_dmg_minus.bed")
    #extend_intervals("accessibility_sorted.bed", "accessibility_extended.bed", 300)

    parser = argparse.ArgumentParser(description="A script that accepts command-line inputs")

    parser.add_argument("tf_name", type=str, help="TF Name")
    parser.add_argument("run_id", type=int, help="Number of Runs")
    args = parser.parse_args()

    start = time.time()

    global genome
    genome = read_fasta("../../data/genome/hg19.fa")

    print(f"Takes {time.time() - start:.2f} seconds to read hg38")

    start = time.time()

    analyze(tf_name=args.tf_name, num_runs=1000, run_id=args.run_id)

    print(f"Takes {time.time() - start:.2f} seconds to process {args.tf_name} - Preprocessed Version")

    
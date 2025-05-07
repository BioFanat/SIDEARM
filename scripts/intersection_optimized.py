import os
import psutil
import math
import time
import struct
import argparse

import numpy as np
import pandas as pd
from intervaltree import IntervalTree, Interval
import statsmodels.stats.multitest as smstats
import matplotlib.pyplot as plt
import logomaker as lm
import scipy.stats as stats
from multiprocessing import Pool
from functools import partial
from pathlib import Path


# =====================================================================================
#                               1. UTILITY / DEBUGGING
# =====================================================================================

def get_memory_usage(label=""):
    """
    Print memory usage of the current process (in MB) for debugging.
    
    Args:
        label (str): Optional string to include in the output.
    """
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 * 1024)
    print(f"[MEMORY] {label} - PID {os.getpid()}: {mem:.2f} MB")


# =====================================================================================
#                               2. FASTA READING
# =====================================================================================

def read_fasta(fasta_file):
    """
    Loads a FASTA file of an entire genome into a dict: {chrom: sequence_str}.
    
    Args:
        fasta_file (str): Path to a reference genome in FASTA format.
    
    Returns:
        dict: Keys = chromosome names, Values = entire chromosome sequence (uppercased).
    """
    get_memory_usage("Before reading FASTA")
    genome_dict = {}
    with open(fasta_file, 'r') as f:
        chrom = None
        seq = []
        for line in f:
            if line.startswith('>'):
                if chrom is not None:
                    genome_dict[chrom] = "".join(seq).upper()
                chrom = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        # Final chromosome
        if chrom is not None:
            genome_dict[chrom] = "".join(seq).upper()
    get_memory_usage("After reading FASTA")
    return genome_dict


# =====================================================================================
#                            3. READING BED/INTERVAL FILES
# =====================================================================================

def read_bed_file(bed_path, mode="damage"):
    """
    Reads a BED-like file. Minimal fields:
    
    If mode='damage', each line is: chrom, start, end, value, strand
    If mode='mutation', each line is: chrom, start [plus optional fields], no explicit strand needed.
    
    Args:
        bed_path (str): Path to bed file
        mode (str): 'damage' or 'mutation'
    
    Returns:
        list of lists: each sub-list is a row of relevant parsed data
    """
    rows = []
    with open(bed_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])

            if mode == "damage":
                # Expect: chrom, start, end, value, strand
                if len(parts) < 5:
                    continue
                end = int(parts[2])
                val = int(parts[3])
                strand = parts[4]
                rows.append([chrom, start, end, val, strand])
            else:
                # mutation mode: at least chrom, start
                # user can decide how to parse the "value" or keep it as 1 for each line.
                # For simplicity, treat each line as one mutation event of value=1
                # (or if you have a 4th column for "value", parse it).
                val = 1
                rows.append([chrom, start, val])
    return rows


def read_tf_bed(tf_bed_file):
    """
    Reads a TF BED file with columns: chrom, start, end, strand (or at least these).
    Returns a list of [chrom, start, end, strand].
    """
    results = []
    with open(tf_bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end   = int(fields[2])
            strand= fields[3]
            results.append([chrom, start, end, strand])
    return results


def hold_interval_indices(interval_file):
    """
    Reads an interval file for decoding simulations, storing (chrom, start) for each line.
    Typically for e.g. /usr/xtmp/bc301/sim_uv_cpd_full/acc_run_full_#.bin
    
    Args:
        interval_file (str): path to intervals bed
    
    Returns:
        list of (chrom, start)
    """
    intervals = []
    with open(interval_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom, start = parts[0], int(parts[1])
            intervals.append((chrom, start))
    return intervals


# =====================================================================================
#                        4. BUILDING THE TF IntervalTree
# =====================================================================================

def build_tf_tree(tf_records, window):
    """
    Given a list of TF motif records [chrom, start, end, strand],
    build an IntervalTree for each chromosome that extends from the motif center ± window.
    
    Args:
        tf_records (list): each entry [chrom, start, end, strand]
        window (int): extension from motif center
    
    Returns:
        dict: {chrom: IntervalTree}, interval data = strand
        int: the first motif length encountered (for plotting).
    """
    trees = {}
    if not tf_records:
        return trees, 0
    # Assume all motifs have the same length. Use the first one as reference:
    motif_len = tf_records[0][2] - tf_records[0][1]

    for rec in tf_records:
        chrom, start, end, strand = rec
        midpoint = (start + end - 1) / 2
        midpoint = math.floor(midpoint) if strand == "+" else math.ceil(midpoint)

        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom].add(Interval(midpoint - window, midpoint + window, strand))
    return trees, motif_len


# =====================================================================================
#                     5. DECODING SIMULATION FILES (DAMAGE vs. MUTATION)
# =====================================================================================

def decode_uint32(code_val, interval_map, mode="damage"):
    """
    Decode a 32-bit integer from a simulation .bin file, returning:
    
    If mode='damage':
      - ('NEW', [some 4-mer code]) or
      - (chrom, start, end, value, strand)
    
    If mode='mutation':
      - ('NEW', [block_id]) or
      - (chrom, pos, value)
    
    We distinguish by the first bit: if it's 1 => 'NEW' marker. If it's 0 => real event.
    
    Args:
        code_val (int): the 32-bit integer
        interval_map (list of (chrom, start)): for region decoding
        mode (str): 'damage' or 'mutation'
    
    Returns:
        tuple: either ('NEW', int/str) or the full decoded record
    """
    bstr = f'{code_val:032b}'
    if bstr[0] == '1':
        # "NEW" marker
        if mode == "damage":
            # last 6 bits store 4-mer info
            fourmer_code = int(bstr[-6:], 2)
            return ("NEW", fourmer_code)
        else:
            # mutation: last 4 bits store block marker
            block_id = int(bstr[-4:], 2)
            return ("NEW", block_id)

    # Otherwise parse a real site
    if mode == "damage":
        # bits: [1:18] => region, [18:27] => pos, bit [27] => strand, [28:32] => value
        region = int(bstr[1:18], 2)
        pos_offset = int(bstr[18:27], 2)
        strand_bit = int(bstr[27], 2)
        dmg_value = int(bstr[28:], 2)

        if region >= len(interval_map):
            return ("NEW", "ERROR")  # out of range
        chrom, region_start = interval_map[region]
        actual_pos = region_start + pos_offset
        strand = '+' if strand_bit == 0 else '-'
        return (chrom, actual_pos, actual_pos+1, dmg_value, strand)

    else:
        # mode='mutation'
        # bits: [1:17] => region, [17:26] => pos_offset, [26:32] => value
        region = int(bstr[1:17], 2)
        pos_offset = int(bstr[17:26], 2)
        mut_val = int(bstr[26:], 2)

        if region >= len(interval_map):
            return ("NEW", -999)
        chrom, region_start = interval_map[region]
        actual_pos = region_start + pos_offset
        return (chrom, actual_pos, mut_val)


# =====================================================================================
#                     6. OVERLAP LOGIC (DAMAGE vs. MUTATION)
# =====================================================================================

def overlap_event(record, tf_tree, genome, mode="damage"):
    """
    For a single *decoded* record, overlap with the TF IntervalTree.
    
    If mode='damage', record = (chrom, start, end, val, strand),
       we gather distance from motif center, plus check same/diff strand.

    If mode='mutation', record = (chrom, pos, val),
       we gather distance from motif center (unstranded).
    
    Args:
        record (tuple): Decoded site from either damage or mutation
        tf_tree (dict): {chrom: IntervalTree}
        genome (dict): reference genome dict (chrom->sequence)
        mode (str): "damage" or "mutation"
    
    Returns:
        list of tuples, each describing an overlap:
        
        - if mode='damage': (distance, value, same_strand_bool)
        - if mode='mutation': (distance, value)
    """
    if mode == "damage":
        chrom, start, end, dmg_val, dmg_strand = record
        if dmg_val <= 0:
            return []
        if chrom not in tf_tree:
            return []
        overlapping = tf_tree[chrom].overlap(start, end)
        results = []
        for iv in overlapping:
            mid = (iv.begin + iv.end) // 2
            tf_strand = iv.data
            # if TF is '-', measure distance from the "end" side
            dist = (end - mid) if (tf_strand == '-') else (start - mid)
            same_strand = (dmg_strand == tf_strand)
            results.append((dist, dmg_val, same_strand))
        return results

    else:
        # mode='mutation'
        chrom, pos, mut_val = record
        if mut_val <= 0:
            return []
        if chrom not in tf_tree:
            return []
        overlapping = tf_tree[chrom].at(pos)
        results = []
        for iv in overlapping:
            mid = (iv.begin + iv.end) // 2
            dist = (pos - mid)
            results.append((dist, mut_val))
        return results


# =====================================================================================
#           7A. PROCESSING A SINGLE SIMULATION RUN  (damage or mutation)
# =====================================================================================

def simulate_run(run_id, interval_map, tf_tree, genome, window, mode="damage"):
    """
    Process a single .bin file for either damage or mutation simulation.
    
    1) Reads the .bin file (acc_run_full_{id}.bin or whichever path logic).
    2) Decodes each 32-bit integer.
    3) For real events, calls overlap_event().
    4) Accumulates counts in a numeric array.
    
    If mode='damage':
       We store 2*(2*window+1) slots: first (2*window+1) for same-strand, next (2*window+1) for opposite-strand.
    If mode='mutation':
       We store (2*window - 1) slots for positions [-window+1..window-1], unstranded.
    
    Args:
        run_id (int): The ID to locate the .bin file
        interval_map (list): from hold_interval_indices
        tf_tree (dict): from build_tf_tree
        genome (dict): reference genome
        window (int): motif extension
        mode (str): "damage" or "mutation"
    
    Returns:
        np.array: aggregated counts.
                 shape = 2*(2*window+1) for damage,
                 shape = (2*window -1) for mutation.
    """
    # Build file path logic
    if mode == "damage":
        file_path = f'/usr/xtmp/bc301/sim_uv_cpd_full/acc_run_full_{run_id}.bin'
        n_slots = 2 * (2 * window + 1)
    else:
        # mode='mutation'
        # example: 'actual' vs 'potential' could also be sub-modes, but we keep it simple here:
        file_path = f'/usr/xtmp/bc301/sim_data_skin_mut_atac/acc_run_{run_id}.bin'
        # or   file_path = f'/usr/xtmp/bc301/sim_data_potential_mut/acc_run_{run_id}.bin'
        n_slots = (2 * window) - 1

    # Prepare array
    accum = np.zeros(n_slots, dtype=float)

    with open(file_path, 'rb') as f:
        raw_data = f.read()
    uint32_arr = np.frombuffer(raw_data, dtype=np.uint32)

    first_new = True
    for code_val in uint32_arr:
        decoded = decode_uint32(code_val, interval_map, mode=mode)
        if decoded[0] == "NEW":
            # indicates new block or end sentinel
            # Usually "NEW" + 0 => end for mutation, "NEW" + '000' => end for damage
            # We'll just check if the block ID is 0 for a sentinel
            if mode == "damage":
                if decoded[1] == 0 and (not first_new):
                    break
            else:
                if decoded[1] == 0 and (not first_new):
                    break
            first_new = False
            continue

        # Real site => check overlap
        overlaps = overlap_event(decoded, tf_tree, genome, mode)
        for ovp in overlaps:
            if mode == "damage":
                dist, val, same_strand = ovp
                idx = dist + window  # shift so index 0 => -window
                if 0 <= idx < (2 * window + 1):
                    if same_strand:
                        # same-strand chunk
                        accum[idx] += val
                    else:
                        # diff-strand chunk => offset by (2*window+1)
                        accum[idx + (2*window + 1)] += val
            else:
                dist, val = ovp
                idx = dist + window - 1  # shift so index 0 => -(window-1)
                if 0 <= idx < n_slots:
                    accum[idx] += val

    return accum


# =====================================================================================
#          7B. PROCESSING ALL SIMULATIONS IN PARALLEL  (damage or mutation)
# =====================================================================================

def build_simulations(num_runs, interval_map, tf_tree, genome, window, mode="damage"):
    """
    Build parallel simulation results for either damage or mutation.

    Args:
        num_runs (int): number of runs
        interval_map (list): from hold_interval_indices
        tf_tree (dict)
        genome (dict)
        window (int)
        mode (str): "damage" or "mutation"
    
    Returns:
        pd.DataFrame: index = positions (float or int), columns = [Run 1..N, min, max, median, mean]
    """
    if mode == "damage":
        # For damage, the final array from each run = 2*(2*window+1)
        # The index will be [-window+0.5..+window+0.5] for the first chunk (same-strand),
        # then the same for diff-strand.
        positions = [i - window + 0.5 for i in range(2*window+1)] * 2
    else:
        # mode='mutation'
        # final array from each run = (2*window - 1)
        # index will be [-window+1..window-1]
        positions = list(range(-window+1, window))

    func = partial(simulate_run,
                   interval_map=interval_map,
                   tf_tree=tf_tree,
                   genome=genome,
                   window=window,
                   mode=mode)

    with Pool() as p:
        all_runs = p.map(func, range(1, num_runs+1))

    # all_runs is a list of np arrays of length n_slots. We want shape = [n_slots, num_runs].
    df = pd.DataFrame(all_runs).transpose()
    df.index = positions
    df.index.name = "Pos"
    df.columns = [f'Run {i}' for i in range(1, num_runs+1)]

    # summary stats
    df['min'] = df.min(axis=1)
    df['max'] = df.max(axis=1)
    df['median'] = df.median(axis=1)
    df['mean'] = df.mean(axis=1)

    return df


# =====================================================================================
#           8A. PROCESSING OBSERVED (REAL) DATA  (damage or mutation)
# =====================================================================================

def process_observed_bed(bed_paths, tf_tree, genome, window, mode="damage"):
    """
    Process one or two BED files for observed data and overlap with TF intervals.

    - If mode='damage', bed_paths should be [plus_path, minus_path].
      We accumulate in a 2*(2*window+1) vector.
    - If mode='mutation', bed_paths should be a single file in the list.
      We accumulate in a (2*window-1) vector.

    Args:
        bed_paths (list): paths to bed files (1 or 2)
        tf_tree (dict): from build_tf_tree
        genome (dict)
        window (int)
        mode (str): "damage" or "mutation"

    Returns:
        pd.DataFrame with columns = (if damage) ["Pos","Dmg","Strand"] or
                                     (if mutation) ["Pos","Mut"].
                 index = "Pos".
    """
    if mode == "damage":
        if len(bed_paths) < 2:
            raise ValueError("For damage mode, need two files: plus_bed, minus_bed.")
        n_slots = 2*(2*window + 1)
        # We store [0..(2*window)] => same strand, then [0..(2*window)] => diff strand
        same_counts = np.zeros(2*window+1, dtype=float)
        diff_counts = np.zeros(2*window+1, dtype=float)

        # Helper to parse and accumulate
        def accumulate_damage(dmg_list):
            for rec in dmg_list:
                chrom, start, end, val, strand = rec
                if val <= 0:
                    continue
                if chrom not in tf_tree:
                    continue
                ovps = overlap_event((chrom, start, end, val, strand),
                                     tf_tree, genome, mode="damage")
                for dist, dval, same_strand in ovps:
                    idx = dist + window
                    if 0 <= idx < (2*window+1):
                        if same_strand:
                            same_counts[idx] += dval
                        else:
                            diff_counts[idx] += dval

        # plus bed
        plus_data = read_bed_file(bed_paths[0], mode="damage")
        accumulate_damage(plus_data)
        # minus bed
        minus_data = read_bed_file(bed_paths[1], mode="damage")
        accumulate_damage(minus_data)

        # Build final DataFrame
        positions = [i - window + 0.5 for i in range(2*window+1)]
        same_rows = [[p, c, "Same"] for p, c in zip(positions, same_counts)]
        diff_rows = [[p, c, "Diff"] for p, c in zip(positions, diff_counts)]
        df = pd.DataFrame(same_rows + diff_rows, columns=["Pos","Dmg","Strand"])
        return df.set_index("Pos")

    else:
        # mode='mutation'
        if len(bed_paths) < 1:
            raise ValueError("For mutation mode, need a single bed file in bed_paths.")
        n_slots = (2*window - 1)
        accum = np.zeros(n_slots, dtype=float)

        # read the bed
        mut_list = read_bed_file(bed_paths[0], mode="mutation")
        for rec in mut_list:
            chrom, pos, val = rec
            if val <= 0:
                continue
            if chrom not in tf_tree:
                continue
            ovps = overlap_event((chrom, pos, val), tf_tree, genome, mode="mutation")
            for dist, dval in ovps:
                idx = dist + window - 1
                if 0 <= idx < n_slots:
                    accum[idx] += dval
        positions = list(range(-window+1, window))
        df = pd.DataFrame({"Pos":positions, "Mut":accum}).set_index("Pos")
        return df


# =====================================================================================
#           8B. SINGLE FUNCTION FOR OBSERVED+SIM MERGE + P-VALUE CALC
# =====================================================================================

def merge_and_add_stats(observed_df, sim_df, motif_len, num_runs, mode="damage"):
    """
    Merge observed data with simulation runs, compute scaling using flanks,
    compute empirical p-values with FDR correction, and store z-scores, etc.

    Observed vs. Sim mapping:
      - if mode='damage': observed_df has columns "Dmg" + "Strand" (optional)
      - if mode='mutation': observed_df has column "Mut".

    Args:
        observed_df (pd.DataFrame): index='Pos'
        sim_df (pd.DataFrame): same index. columns = runs + [min, max, median, mean].
        motif_len (int): length of the motif (used for defining core region).
        num_runs (int)
        mode (str): 'damage' or 'mutation'

    Returns:
        pd.DataFrame: includes the observed col (Dmg or Mut), the simulation runs,
                      plus 'P-Value-Top','P-Value-Bottom','z', etc.
                      Rows are re-indexed by 'Pos' (potentially for each strand if damage).
    """
    df = pd.concat([observed_df, sim_df], axis=1)
    df.index.name = "Pos"
    df.reset_index(inplace=True)

    # Which col do we use? "Dmg" or "Mut"?
    obs_col = "Dmg" if mode == "damage" else "Mut"

    # We'll define the motif range for possible scaling using flanks
    half_len_left = -(math.ceil(motif_len / 2) - 1)
    half_len_right = motif_len // 2

    # Define flank masks
    if mode == "damage":
        # Because for damage we have 'Strand' in df, we do separate scaling for each strand half.
        return _add_stats_damage(df, obs_col, num_runs, half_len_left, half_len_right)
    else:
        # mode='mutation'
        return _add_stats_mutation(df, obs_col, num_runs, half_len_left, half_len_right)


def _add_stats_damage(df, obs_col, num_runs, left_motif_edge, right_motif_edge):
    """
    Internal subroutine: damage p-value calc with same/diff strand scaling.
    Each half (left vs. right) is scaled separately *per strand* using flanks.
    """
    run_cols = slice(3, -4)  # after index=Pos, obs_col=1, Strand=2 => runs up to the last 4 stats
    df_runs = df.iloc[:, run_cols]
    # For each Strand, do left- and right-flank scaling:
    def compute_factor(subset):
        # ratio = (mean observed) / (mean of sim means)
        return (subset[obs_col].mean()) / (subset['mean'].mean())

    # Mark flanks
    df['LeftFlank'] = (df['Pos'] < (left_motif_edge - 5))
    df['RightFlank'] = (df['Pos'] > (right_motif_edge + 5))
    df['LeftHalf'] = (df['Pos'] < 0)
    df['RightHalf'] = (df['Pos'] >= 0)

    # Collect scaling factors in dict: e.g. factors[(strand, 'left')] = 1.2
    factors = {}
    for strand_val in df['Strand'].unique():
        strand_mask = (df['Strand'] == strand_val)
        left_flank_data = df[strand_mask & df['LeftFlank']]
        right_flank_data = df[strand_mask & df['RightFlank']]

        if not left_flank_data.empty:
            factors[(strand_val, 'left')] = compute_factor(left_flank_data)
        if not right_flank_data.empty:
            factors[(strand_val, 'right')] = compute_factor(right_flank_data)

    # Apply scaling
    for strand_val in df['Strand'].unique():
        # left side
        mask_left = (df['Strand'] == strand_val) & df['LeftHalf']
        if (strand_val, 'left') in factors:
            factor_l = factors[(strand_val, 'left')]
            df.loc[mask_left, df.columns[run_cols]] *= factor_l
        # right side
        mask_right = (df['Strand'] == strand_val) & df['RightHalf']
        if (strand_val, 'right') in factors:
            factor_r = factors[(strand_val, 'right')]
            df.loc[mask_right, df.columns[run_cols]] *= factor_r

    # Recompute summary stats
    df_runs = df.iloc[:, run_cols]  # after scaling
    df['min'] = df_runs.min(axis=1)
    df['max'] = df_runs.max(axis=1)
    df['median'] = df_runs.median(axis=1)
    df['mean'] = df_runs.mean(axis=1)
    df['std'] = df_runs.std(axis=1)

    # Compute empirical p-values
    pvals_top = []
    pvals_bottom = []
    obs_vals = df[obs_col].values
    for i in range(len(df)):
        sim_vals = df_runs.iloc[i].values
        obs = obs_vals[i]
        if np.allclose(sim_vals, obs):
            pvals_top.append(1.0)
            pvals_bottom.append(1.0)
        else:
            pvals_top.append((sim_vals > obs).mean())
            pvals_bottom.append((sim_vals < obs).mean())

    q_top = smstats.fdrcorrection(pvals_top)[1]
    q_bottom = smstats.fdrcorrection(pvals_bottom)[1]
    df['P-Value-Top'] = q_top
    df['P-Value-Bottom'] = q_bottom

    # z-scores
    df['z'] = (df[obs_col] - df['mean']) / (df['std'] + 1e-9)

    df.drop(columns=['LeftFlank','RightFlank','LeftHalf','RightHalf'], inplace=True)
    df.set_index("Pos", inplace=True)
    return df


def _add_stats_mutation(df, obs_col, num_runs, left_motif_edge, right_motif_edge):
    """
    Internal subroutine: unstranded 'mutation' p-value calc. Single half-based scaling.
    """
    run_cols = slice(2, -4)  # after index=Pos, obs_col => runs up to last 4 stats
    df_runs = df.iloc[:, run_cols]

    # Mark flanks
    df['LeftFlank'] = (df['Pos'] < (left_motif_edge - 5))
    df['RightFlank'] = (df['Pos'] > (right_motif_edge + 5))
    df['LeftHalf'] = (df['Pos'] < 0)
    df['RightHalf'] = (df['Pos'] >= 0)

    # Compute factors
    def compute_factor(subdf):
        return (subdf[obs_col].mean()) / (subdf['mean'].mean())

    left_data = df[df['LeftFlank']]
    right_data = df[df['RightFlank']]
    factors = {}
    if not left_data.empty:
        factors['left'] = compute_factor(left_data)
    if not right_data.empty:
        factors['right'] = compute_factor(right_data)

    # Apply
    mask_left = df['LeftHalf']
    if 'left' in factors:
        df.loc[mask_left, df.columns[run_cols]] *= factors['left']
    mask_right = df['RightHalf']
    if 'right' in factors:
        df.loc[mask_right, df.columns[run_cols]] *= factors['right']

    # Recompute stats
    df_runs = df.iloc[:, run_cols]
    df['min'] = df_runs.min(axis=1)
    df['max'] = df_runs.max(axis=1)
    df['median'] = df_runs.median(axis=1)
    df['mean'] = df_runs.mean(axis=1)
    df['std'] = df_runs.std(axis=1)

    # P-value
    pvals_top = []
    pvals_bottom = []
    obs_vals = df[obs_col].values
    for i in range(len(df)):
        sim_vals = df_runs.iloc[i].values
        obs = obs_vals[i]
        if np.allclose(sim_vals, obs):
            pvals_top.append(1.0)
            pvals_bottom.append(1.0)
        else:
            pvals_top.append((sim_vals > obs).mean())
            pvals_bottom.append((sim_vals < obs).mean())

    q_top = smstats.fdrcorrection(pvals_top)[1]
    q_bottom = smstats.fdrcorrection(pvals_bottom)[1]
    df['P-Value-Top'] = q_top
    df['P-Value-Bottom'] = q_bottom

    df['z'] = (df[obs_col] - df['mean']) / (df['std'] + 1e-9)
    df.drop(columns=['LeftFlank','RightFlank','LeftHalf','RightHalf'], inplace=True)
    df.set_index("Pos", inplace=True)
    return df


# =====================================================================================
#                           9. LOGO & PLOTTING
# =====================================================================================

def rev_complement(seq):
    """Simple reverse complement of a DNA sequence."""
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join(comp.get(b,b) for b in seq[::-1])


def build_logo_sequences(observed_df, bed_paths, tf_tree, genome, window, mode="damage"):
    """
    Optionally extract the 31bp sequences around each TF center for real data (observed).
    If desired, one can do the same for simulated data, but that often is huge.

    For each overlap, we store the 31mer around the motif center.

    Returns:
        plus_seq_list, minus_seq_list
    """
    # Typically you only do big sequence extraction for the actual data, not for all 10k sims.
    plus_list = []
    minus_list = []

    if mode == "damage":
        # We'll re-run the overlaps or we can incorporate in the code above for speed.
        # For brevity, we just show a placeholder approach.
        pass
    else:
        # same placeholder for mutation
        pass
    
    # Implementation left as an exercise, but typically you'd parse the bed events,
    # check tf_tree overlap, find each center => genome[center -15 : center +16], etc.
    # Then store in plus_list / minus_list according to the motif strand.
    return plus_list, minus_list


def generate_differential_logo(plus_seqs, minus_seqs, motif_len, ax):
    """
    Builds a difference sequence logo (plus minus).
    Expects plus_seqs and minus_seqs all of length 2*window+1, e.g. 31.

    Args:
        plus_seqs (list of str)
        minus_seqs (list of str)
        motif_len (int)
        ax (matplotlib.Axes)
    """
    import logomaker as lm

    # Filter out sequences that do not have the correct length
    # ...
    # Then build plus and minus logos -> subtract -> plot
    pass  # Implementation similar to prior snippet


def plot_stacked(df, tf_name, motif_len, window, mode="damage", output_prefix="output"):
    """
    Creates a multi-panel figure showing:
      1) Observed vs. sim envelope (+ significance markers)
      2) -log10(q-values)
      3) Z-score
      4) Sequence logo (optional)
    
    Args:
        df (pd.DataFrame): The final merged data with columns ["Dmg" or "Mut", "P-Value-Top", "P-Value-Bottom", "z", etc.]
        tf_name (str)
        motif_len (int)
        window (int)
        mode (str): "damage" or "mutation"
        output_prefix (str)
    """
    # Implementation is analogous to the prior stacked plotting logic.
    # Consolidate for both modes, e.g. if "damage" => color by strand; if "mutation" => single color.
    # Then produce a final .png
    pass


# =====================================================================================
#           10. WRITING SIGNIFICANCE RESULTS  (damage or mutation)
# =====================================================================================

def write_significance_file(df, tf_name, mode="damage", alpha=0.05, out_dir="output"):
    """
    Writes a text file summarizing significance calls (positions in [-15..15], etc.).
    For damage: includes Strand; for mutation: no Strand.
    
    Args:
        df (pd.DataFrame): final results. Columns = "P-Value-Top", "P-Value-Bottom", "z", "Strand"? ...
        tf_name (str)
        mode (str)
        alpha (float)
        out_dir (str)
    """
    out_path = Path(out_dir) / f"{tf_name}_{mode}.out"
    # Filter to -15..15
    subset = df[(df.index >= -15) & (df.index <= 15)]

    if mode == "damage":
        # We have columns: "Dmg", "Strand", "P-Value-Top", "P-Value-Bottom", "z"
        min_pval_top = subset["P-Value-Top"].min()
        min_pval_bottom = subset["P-Value-Bottom"].min()
        sig_top = subset[subset["P-Value-Top"] < alpha]
        sig_bot = subset[subset["P-Value-Bottom"] < alpha]

        zmax = subset["z"].max()
        zmin = subset["z"].min()

        with open(out_path, 'w') as f:
            f.write(f"SIGNIFICANCE ANALYSIS FOR {tf_name} ({mode.upper()})\n\n")
            f.write("MinP-Top\tMinP-Bottom\tMaxZ\tMinZ\n")
            f.write(f"{min_pval_top}\t{min_pval_bottom}\t{zmax}\t{zmin}\n\n")
            f.write("Significant (Enriched):\n")
            for pos, row in sig_top.iterrows():
                f.write(f"{pos}\t{row['Strand']}\t{row['P-Value-Top']}\t{row['z']}\n")
            f.write("\nSignificant (Protected):\n")
            for pos, row in sig_bot.iterrows():
                f.write(f"{pos}\t{row['Strand']}\t{row['P-Value-Bottom']}\t{row['z']}\n")
    else:
        # mode='mutation'
        # columns: "Mut", "P-Value-Top", "P-Value-Bottom", "z"
        min_pval_top = subset["P-Value-Top"].min()
        min_pval_bottom = subset["P-Value-Bottom"].min()
        sig_top = subset[subset["P-Value-Top"] < alpha]
        sig_bot = subset[subset["P-Value-Bottom"] < alpha]
        zmax = subset["z"].max()
        zmin = subset["z"].min()

        with open(out_path, 'w') as f:
            f.write(f"SIGNIFICANCE ANALYSIS FOR {tf_name} ({mode.upper()})\n\n")
            f.write("MinP-Top\tMinP-Bottom\tMaxZ\tMinZ\n")
            f.write(f"{min_pval_top}\t{min_pval_bottom}\t{zmax}\t{zmin}\n\n")
            f.write("Significant (Increased Mutations):\n")
            for pos, row in sig_top.iterrows():
                f.write(f"{pos}\t{row['P-Value-Top']}\t{row['z']}\n")
            f.write("\nSignificant (Decreased Mutations):\n")
            for pos, row in sig_bot.iterrows():
                f.write(f"{pos}\t{row['P-Value-Bottom']}\t{row['z']}\n")


# =====================================================================================
#                                11. MAIN PIPELINE
# =====================================================================================

def analyze_tf(
    tf_name, 
    mode, 
    genome_fasta, 
    tf_plus_bed, 
    tf_minus_bed, 
    interval_bed, 
    observed_beds, 
    num_runs=1000, 
    window=15,
    output_dir="output"
):
    """
    High-level pipeline that:
      1) Loads genome.
      2) Reads TF sites (plus + minus), build IntervalTree.
      3) Processes observed data (plus/minus if damage, single bed if mutation).
      4) Runs parallel simulation for num_runs.
      5) Merge, scale, p-value, etc.
      6) Write .out significance file + optional figure.
    
    Args:
        tf_name (str)
        mode (str): 'damage' or 'mutation'
        genome_fasta (str): path to reference genome
        tf_plus_bed (str): bed for TF plus
        tf_minus_bed (str): bed for TF minus
        interval_bed (str): bed for intervals (region map for decoding)
        observed_beds (list): for damage => [plus_path, minus_path], for mutation => [single_path]
        num_runs (int)
        window (int)
        output_dir (str)
    """
    # 1) Load genome
    genome = read_fasta(genome_fasta)

    # 2) Build TF records => IntervalTree
    plus_tf = read_tf_bed(tf_plus_bed)
    minus_tf= read_tf_bed(tf_minus_bed)
    tf_records = plus_tf + minus_tf
    tf_tree, motif_len = build_tf_tree(tf_records, window)

    # 3) Process observed data
    obs_df = process_observed_bed(observed_beds, tf_tree, genome, window, mode=mode)

    # 4) Build simulations
    interval_map = hold_interval_indices(interval_bed)
    sim_df = build_simulations(num_runs, interval_map, tf_tree, genome, window, mode=mode)

    # 5) Merge + Stats
    final_df = merge_and_add_stats(obs_df, sim_df, motif_len, num_runs, mode=mode)

    # 6) Write significance
    out_csv_path = Path(output_dir) / f"{tf_name}_{mode}_final.csv"
    final_df.to_csv(out_csv_path)
    write_significance_file(final_df, tf_name, mode=mode, out_dir=output_dir)

    # 7) (Optional) Make figure
    # plot_stacked(final_df, tf_name, motif_len, window, mode=mode, output_prefix=output_dir)
    # or create your own custom multi-panel figure.

    print(f"Done analyzing {tf_name} ({mode}); results saved to {out_csv_path}")


# =====================================================================================
#                               12. CLI ENTRY POINT
# =====================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Refactored UV CPD / Mutation Analysis")
    parser.add_argument("tf_name", type=str, help="TF name")
    parser.add_argument("mode", choices=["damage","mutation"], help="Analysis mode")
    parser.add_argument("runs", type=int, help="Number of simulation runs")
    parser.add_argument("--window", type=int, default=20, help="Window around motif center")
    args = parser.parse_args()

    # Example usage below: you would update to match your actual file paths.
    # Or you can pass them in from the command line or a config file.
    # For demonstration, we just show placeholders:

    genome_fasta = "../../data/genome/hg19.fa"
    tf_plus_bed  = f"../../data/tfbs/fib_atac/{args.tf_name}/{args.tf_name}_plus_high.bed"
    tf_minus_bed = f"../../data/tfbs/fib_atac/{args.tf_name}/{args.tf_name}_minus_high.bed"
    interval_bed = "../../data/raw/atac_150bp_intervals_merged.bed"

    if args.mode == "damage":
        observed_beds = ["../../data/damages/atac_nb_plus.bed", "../../data/damages/atac_nb_minus.bed"]
    else:
        observed_beds = ["../../data/raw/atac_mutations_transitions_C_only.bed"]

    start_time = time.time()
    analyze_tf(
        tf_name=args.tf_name,
        mode=args.mode,
        genome_fasta=genome_fasta,
        tf_plus_bed=tf_plus_bed,
        tf_minus_bed=tf_minus_bed,
        interval_bed=interval_bed,
        observed_beds=observed_beds,
        num_runs=args.runs,
        window=args.window,
        output_dir="../../output/messing_around"
    )
    print(f"Total time: {time.time() - start_time:.2f} s")
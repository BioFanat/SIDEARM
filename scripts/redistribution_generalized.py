#!/usr/bin/env python3
"""
Generalized redistribution for SIDEARM

This script redistributes damage counts across k-mers using a configurable
framework for different damage types.
"""

import argparse
import os
import struct
import time
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import pandas as pd

# Import our new configuration system
from kmer_config import KmerConfig, get_config_by_name, get_default_config

###############################################################################
#                       1) PARSING THE INPUT .BIN FILE                        #
###############################################################################
def parse_input_bin(file_path: str, config: KmerConfig) -> Dict[int, np.ndarray]:
    """
    Parse an input binary file into a dictionary of arrays.
    
    Args:
        file_path: Path to the input binary file
        config: KmerConfig object for decoding
        
    Returns:
        Dictionary mapping ID to numpy array of (region, pos, strand_bit) tuples
    """
    data_dict = {}
    with open(file_path, 'rb') as bin_file:
        current_id = None
        block_rows = []  # temporary list for building the array

        while True:
            chunk = bin_file.read(4)
            if len(chunk) < 4:
                # End of file => store the last block if any
                if current_id is not None and block_rows:
                    arr_np = np.array(block_rows, dtype=np.uint16)
                    data_dict[current_id] = arr_np
                break

            val = struct.unpack('I', chunk)[0]
            decoded = config.decode_site(val)

            if decoded[0] == "NEW":
                # Start of a new ID block
                # => first store the old block, if it exists
                if current_id is not None and block_rows:
                    arr_np = np.array(block_rows, dtype=np.uint16)
                    data_dict[current_id] = arr_np

                # Reset for the new ID
                current_id = decoded[1]
                block_rows = []

            else:
                # Normal line => decode fields (region, pos, strand, dmg)
                site_type, region_val, pos_val, strand_sym, _ = decoded
                strand_bit = 0 if strand_sym == '+' else 1
                block_rows.append((region_val, pos_val, strand_bit))

        # (End while) we already stored the last block above if any
    return data_dict

###############################################################################
#                          2) DAMAGE REDISTRIBUTION                           #
###############################################################################
def redistribute_uniform(array_3col: np.ndarray, total_dmg: int) -> np.ndarray:
    """
    Distribute damage uniformly across sites.
    
    Args:
        array_3col: Numpy array of (region, pos, strand_bit) tuples
        total_dmg: Total damage to distribute
        
    Returns:
        Numpy array of (region, pos, strand_bit, damage) tuples
    """
    N = array_3col.shape[0]
    if N == 0 or total_dmg < 1:
        # No sites or no damage => empty
        return np.zeros((0, 4), dtype=np.uint16)

    if N == 1:
        # Single site => all damage there
        return np.array([
            [array_3col[0,0], array_3col[0,1], array_3col[0,2], total_dmg]
        ], dtype=np.uint16)

    # Uniform probabilities
    p = np.ones(N, dtype=float) / N
    counts = np.random.multinomial(total_dmg, p)

    # Filter zero-damage
    nz_mask = (counts > 0)
    if not np.any(nz_mask):
        # It's possible (though unlikely) that the multinomial 
        # gave zero to every site if total_dmg is small. 
        return np.zeros((0, 4), dtype=np.uint16)

    # Build result array
    sub_3col = array_3col[nz_mask]
    dmg_vals = counts[nz_mask].astype(np.uint16)
    result = np.column_stack([sub_3col, dmg_vals])  # shape (M,4)
    return result

###############################################################################
#          3) PERFORM SIMULATIONS (generating multiple run .bin files)        #
###############################################################################
def perform_simulations(
    file_path: str,
    index_path: str,
    runs: int,
    start_run_id: int,
    output_prefix: str,
    damage_type: str = "CPD"
):
    """
    Perform multiple simulation runs.
    
    Args:
        file_path: Path to the input binary file
        index_path: Path to the damage index file
        runs: Number of simulation runs
        start_run_id: Starting run ID
        output_prefix: Prefix for output files
        damage_type: Damage type (default: 'CPD')
    """
    # Get configuration for this damage type
    try:
        config = get_config_by_name(damage_type)
    except ValueError:
        print(f"Warning: Unknown damage type '{damage_type}'. Using default (CPD).")
        config = get_default_config()
    
    print(f"Using configuration for {config.name} with kmer length {config.kmer_length}")

    if runs < 1:
        print("[WARNING] 'runs' must be >=1. Doing nothing.")
        return

    # 1) Parse input bin
    print("[INFO] Parsing input .bin file...")
    data_dict = parse_input_bin(file_path, config)
    all_ids = sorted(data_dict.keys())
    if not all_ids:
        print("[ERROR] No IDs found in input file. Exiting.")
        return
    max_id = max(all_ids)

    # 2) Read index => total damage for each ID
    print("[INFO] Reading index file for total damage info...")
    df = pd.read_csv(index_path, sep="\t", header=None,
                     names=["dmg_plus", "dmg_minus", "counts_plus", "counts_minus"])
    df['dmg_total'] = df['dmg_plus'] + df['dmg_minus']
    dmg_array = df['dmg_total'].values  # 1D array of length = #rows in the file

    # 3) Loop for runs
    print(f"[INFO] Starting {runs} run(s).")
    start_t = time.time()
    np.random.seed(None)

    for i in range(runs):
        # The actual run ID on disk
        run_id = start_run_id*runs + i + 1
        out_fname = os.path.join(output_prefix, f"acc_run_full_{run_id}.bin")

        with open(out_fname, 'wb') as out_f:
            for cur_id in all_ids:
                # Create sentinel for this ID
                sentinel_code = config.encode_sentinel(cur_id)
                out_f.write(struct.pack('I', sentinel_code))

                # Retrieve total damage if the ID < len(dmg_array)
                if 0 <= cur_id < len(dmg_array):
                    total_dmg = int(dmg_array[cur_id])
                else:
                    total_dmg = 0

                # Redistribute
                arr_3col = data_dict[cur_id]
                res_4col = redistribute_uniform(arr_3col, total_dmg)
                # res_4col => shape (M,4): (region, pos, strand_bit, damage)

                # Write
                for (r_val, p_val, strand_bit, dmg_val) in res_4col:
                    # Convert strand_bit => strand symbol
                    strand_char = '+' if strand_bit == 0 else '-'
                    encoded = config.encode_site(r_val, p_val, strand_char, dmg_val)
                    out_f.write(struct.pack('I', encoded))

            # Final sentinel => 'NEW' with (max_id + 1)
            final_sentinel = config.encode_sentinel(max_id + 1)
            out_f.write(struct.pack('I', final_sentinel))

        print(f"   [INFO] Finished run {run_id} => {out_fname}")

    elapsed = time.time() - start_t
    print(f"[INFO] Completed {runs} run(s) in {elapsed:.2f} sec total.")

###############################################################################
#                           4) OPTIONAL DECODER CHECK                          #
###############################################################################
def simple_decode_check(bin_path: str, index_path: str, damage_type: str = "CPD"):
    """
    Check that a binary file matches the expected damage distribution.
    
    Args:
        bin_path: Path to the binary file to check
        index_path: Path to the damage index file
        damage_type: Damage type (default: 'CPD')
    """
    # Get configuration for this damage type
    try:
        config = get_config_by_name(damage_type)
    except ValueError:
        print(f"Warning: Unknown damage type '{damage_type}'. Using default (CPD).")
        config = get_default_config()
    
    print(f"Using configuration for {config.name} with kmer length {config.kmer_length}")
    
    df = pd.read_csv(index_path, sep="\t", header=None,
                     names=["dmg_plus", "dmg_minus", "counts_plus", "counts_minus"])
    df['dmg_total'] = df['dmg_plus'] + df['dmg_minus']
    dmg_array = df['dmg_total'].values

    with open(bin_path, 'rb') as f:
        current_id = None
        total_dmg = 0

        while True:
            chunk = f.read(4)
            if len(chunk) < 4:
                # end
                if current_id is not None:
                    print(f"End => ID {current_id} total damage found = {total_dmg}")
                break
            val = struct.unpack('I', chunk)[0]
            dec = config.decode_site(val)
            if dec[0] == "NEW":
                # If we had an ID already, let's see how it matched
                if current_id is not None:
                    if current_id < len(dmg_array):
                        expected = dmg_array[current_id]
                        msg = f"ID={current_id}: total_dmg={total_dmg}, expected={expected}"
                    else:
                        msg = f"ID={current_id}: total_dmg={total_dmg}, expected=N/A (out of index)"
                    print(msg)

                current_id = dec[1]
                total_dmg = 0
            else:
                # normal site => add damage
                dmg_val = dec[4]  # damage is at position 4 in the tuple
                total_dmg += dmg_val

###############################################################################
#                           5) MAIN (CLI INTERFACE)                           #
###############################################################################
def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Generalized distribution of damages across k-mers.")
    parser.add_argument("--id", type=int, required=True,
                        help="Starting run ID offset (e.g., 1 => first run is '1', second run is '2', etc.)")
    parser.add_argument("--runs", type=int, default=10,
                        help="Number of runs to perform in one go.")
    parser.add_argument("--input", type=str, required=True,
                        help="Path to the input .bin file (k-mer data).")
    parser.add_argument("--index", type=str, required=True,
                        help="Path to the .out or .csv with columns for damage (dmg_total).")
    parser.add_argument("--output", type=str, default=".",
                        help="Directory or prefix where output .bin files will go.")
    parser.add_argument("--damage-type", type=str, default="CPD",
                        help="Damage type (CPD, BPDE, or custom).")
    parser.add_argument("--decode", type=str, default="",
                        help="If provided, decode this file (instead of performing simulations).")

    args = parser.parse_args()

    if args.decode:
        # Just decode-check
        if not os.path.isfile(args.decode):
            print(f"[ERROR] decode file not found: {args.decode}")
        else:
            simple_decode_check(args.decode, args.index, args.damage_type)
    else:
        # Perform simulations
        perform_simulations(
            file_path=args.input,
            index_path=args.index,
            runs=args.runs,
            start_run_id=args.id - 1,
            output_prefix=args.output,
            damage_type=args.damage_type
        )

if __name__ == "__main__":
    main()
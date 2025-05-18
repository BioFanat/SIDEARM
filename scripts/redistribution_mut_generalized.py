#!/usr/bin/env python3
"""
Generalized mutation redistribution for SIDEARM

This script redistributes mutation counts across k-mers using a configurable
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
from kmer_config import KmerConfig, get_config_by_name, get_default_config, BPDE_CONFIG

###############################################################################
#                       1) PARSING THE INPUT .BIN FILE                        #
###############################################################################
def parse_input_bin(bin_path: str, config: KmerConfig) -> Dict[int, np.ndarray]:
    """
    Parse an input binary file into a dictionary of arrays.
    
    Args:
        bin_path: Path to the input binary file
        config: KmerConfig object for decoding
        
    Returns:
        Dictionary mapping ID to numpy array of (region, pos) tuples
    """
    data_dict = {}
    current_id = None
    block_rows = []

    with open(bin_path, 'rb') as f:
        while True:
            chunk = f.read(4)
            if len(chunk) < 4:
                # End of file
                # store the last block
                if current_id is not None and len(block_rows) > 0:
                    arr_np = np.array(block_rows, dtype=np.uint32)
                    data_dict[current_id] = arr_np
                break

            val = struct.unpack('I', chunk)[0]
            decoded = config.decode_site(val)

            if decoded[0] == "NEW":
                # store the previous block
                if current_id is not None and len(block_rows) > 0:
                    arr_np = np.array(block_rows, dtype=np.uint32)
                    data_dict[current_id] = arr_np

                new_id = decoded[1]
                if new_id == 16:  # Typical stopping condition
                    break

                current_id = new_id
                block_rows = []
            else:
                # normal line => (type, region, pos, strand, dmg)
                # but we ignore 'dmg' and 'strand' in the input for the re-distribution step
                _, region_val, pos_val, _, _ = decoded
                block_rows.append([region_val, pos_val])
        
        if current_id is not None and len(block_rows) > 0:
            arr_np = np.array(block_rows, dtype=np.uint32)
            data_dict[current_id] = arr_np

    return data_dict


###############################################################################
#                            2) READ INDEX FILE                               #
###############################################################################
def read_index_file(index_path: str) -> np.ndarray:
    """
    Read mutation counts from an index file.
    
    Args:
        index_path: Path to the index file
        
    Returns:
        Numpy array of mutation counts
    """
    df = pd.read_csv(index_path, sep='\t', header=None)
    # If the first column is the total mutation count:
    total_col = df[0].values  # <— Adjust if necessary

    return total_col


###############################################################################
#                    3) UNIFORM REDISTRIBUTION OF MUTATIONS                   #
###############################################################################
def redistribute_uniform(array_2col: np.ndarray, total_mut: int) -> np.ndarray:
    """
    Distribute mutations uniformly across sites.
    
    Args:
        array_2col: Numpy array of (region, pos) tuples
        total_mut: Total mutations to distribute
        
    Returns:
        Numpy array of (region, pos, damage) tuples
    """
    N = len(array_2col)
    if N == 0 or total_mut < 1:
        return np.zeros((0, 3), dtype=np.uint32)

    if N == 1:
        # all mutations to the single site
        return np.array([[array_2col[0,0], array_2col[0,1], total_mut]], dtype=np.uint32)

    # uniform probabilities = 1/N
    p = np.ones(N) / N
    counts = np.random.multinomial(total_mut, p)
    nz_mask = (counts > 0)

    if not np.any(nz_mask):
        # edge case if total_mut is small and didn't get allocated
        return np.zeros((0, 3), dtype=np.uint32)

    out = np.column_stack([array_2col[nz_mask], counts[nz_mask]])
    return out


###############################################################################
#           4) PERFORM SIMULATIONS, WRITING MULTIPLE .BIN OUTPUTS            #
###############################################################################
def perform_simulations(
    input_bin: str,
    index_path: str,
    runs: int,
    start_run_id: int,
    output_prefix: str,
    damage_type: str = "BPDE"
):
    """
    Perform multiple simulation runs.
    
    Args:
        input_bin: Path to the input binary file
        index_path: Path to the mutation index file
        runs: Number of simulation runs
        start_run_id: Starting run ID
        output_prefix: Prefix for output files
        damage_type: Damage type (default: 'BPDE')
    """
    # Get configuration for this damage type
    try:
        config = get_config_by_name(damage_type)
    except ValueError:
        print(f"Warning: Unknown damage type '{damage_type}'. Using default (BPDE).")
        config = BPDE_CONFIG  # Use BPDE as default for mutations
    
    print(f"Using configuration for {config.name} with kmer length {config.kmer_length}")

    if runs < 1:
        print("[WARNING] 'runs' must be >= 1. No simulations performed.")
        return

    # 1) parse the .bin
    print("[INFO] Parsing input bin...")
    data_dict = parse_input_bin(input_bin, config)
    all_ids = sorted(data_dict.keys())
    if not all_ids:
        print("[ERROR] No IDs found in the input bin. Exiting.")
        return
    max_id = max(all_ids)

    # 2) read the index file
    print("[INFO] Reading index file for total mutations...")
    total_mut_array = read_index_file(index_path)

    # 3) do runs
    print(f"[INFO] Starting {runs} run(s).")
    start_t = time.time()
    np.random.seed(None)  # or set a seed if reproducibility is needed

    for i in range(runs):
        # the actual run ID on disk
        run_id = (start_run_id-1) * runs + (i + 1)
        out_fname = os.path.join(output_prefix, f"mut_run_full_{run_id}.bin")

        with open(out_fname, 'wb') as out_f:
            for cur_id in all_ids:
                # Create sentinel for this ID
                sentinel_code = config.encode_sentinel(cur_id)
                out_f.write(struct.pack('I', sentinel_code))

                # total mutations for cur_id
                if 0 <= cur_id < len(total_mut_array):
                    total_mut = int(total_mut_array[cur_id])
                else:
                    total_mut = 0

                arr_2col = data_dict[cur_id]
                res_3col = redistribute_uniform(arr_2col, total_mut)
                # shape (M,3): (region, pos, damage)

                for row in res_3col:
                    region_val = int(row[0])
                    pos_val = int(row[1])
                    dmg_val = int(row[2])
                    
                    # For mutations, we typically don't use strand
                    encoded = config.encode_site(region_val, pos_val, '+', dmg_val)
                    out_f.write(struct.pack('I', encoded))

            # final sentinel => NEW with (max_id+1) or 0 if you want to end
            final_sentinel = config.encode_sentinel(max_id+1)
            out_f.write(struct.pack('I', final_sentinel))

        print(f"   [INFO] Finished run {i+1} => {out_fname}")

    elapsed = time.time() - start_t
    print(f"[INFO] Completed {runs} run(s) in {elapsed:.2f} sec total.")


###############################################################################
#                            5) OPTIONAL DECODER                              #
###############################################################################
def decode_result(bin_path: str, index_path: str, damage_type: str = "BPDE"):
    """
    A minimal check. We read the generated .bin file and sum the total damage 
    for each ID (only if ID < len(index_array)). Then compare to what's in the index.
    
    Args:
        bin_path: Path to the binary file to check
        index_path: Path to the mutation index file
        damage_type: Damage type (default: 'BPDE')
    """
    # Get configuration for this damage type
    try:
        config = get_config_by_name(damage_type)
    except ValueError:
        print(f"Warning: Unknown damage type '{damage_type}'. Using default (BPDE).")
        config = BPDE_CONFIG  # Use BPDE as default for mutations
    
    print(f"Using configuration for {config.name} with kmer length {config.kmer_length}")
    
    index_array = read_index_file(index_path)

    with open(bin_path, 'rb') as f:
        current_id = None
        total_dmg_for_id = 0

        while True:
            chunk = f.read(4)
            if len(chunk) < 4:
                # end
                if current_id is not None:
                    # final check
                    if current_id < len(index_array):
                        print(f"ID={current_id}: Found total={total_dmg_for_id}, expected={index_array[current_id]}")
                break

            val = struct.unpack('I', chunk)[0]
            dec = config.decode_site(val)
            if dec[0] == "NEW":
                # new block
                # check the old block if existed
                if current_id is not None:
                    if current_id < len(index_array):
                        print(f"ID={current_id}: Found total={total_dmg_for_id}, expected={index_array[current_id]}")

                current_id = dec[1]
                total_dmg_for_id = 0
                if current_id == 16:
                    # typical end sentinel
                    break
            else:
                # normal line => add damage (position 4 in the tuple)
                dmg_val = dec[4]
                total_dmg_for_id += dmg_val


###############################################################################
#                               6) MAIN CLI                                   #
###############################################################################
def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Generalized distribution of mutations across k-mers.")
    parser.add_argument("--input", type=str, required=False,
                        help="Path to the input .bin file with 'NEW' + (region,pos,dmg).")
    parser.add_argument("--index", type=str, required=True,
                        help="Path to the tab-delimited index file containing total mutations per ID row.")
    parser.add_argument("--output", type=str, default=".",
                        help="Output directory (prefix) for .bin results.")
    parser.add_argument("--runs", type=int, default=1,
                        help="Number of simulation runs.")
    parser.add_argument("--start_id", type=int, default=1,
                        help="Numeric ID offset for naming the output runs.")
    parser.add_argument("--damage-type", type=str, default="BPDE",
                        help="Damage type (BPDE or custom).")
    parser.add_argument("--check", type=str, default="",
                        help="If provided, will run a decode check on that .bin instead of performing simulations.")

    args = parser.parse_args()

    if args.check:
        decode_result(args.check, args.index, args.damage_type)
    else:
        if not args.input:
            parser.error("--input is required when not using --check")
            
        perform_simulations(
            input_bin=args.input,
            index_path=args.index,
            runs=args.runs,
            start_run_id=args.start_id,
            output_prefix=args.output,
            damage_type=args.damage_type
        )

if __name__ == "__main__":
    main()
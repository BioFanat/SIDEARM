import argparse
import os
import struct
import time

import numpy as np
import pandas as pd

###############################################################################
#                           1) DECODING / ENCODING                            #
###############################################################################

def decode_nochrom(site_code: int):
    """
    This is exactly the same logic you provided in the original code.

    Decodes a 32-bit integer (site_code) into:
       (region, position, strand, damage)
    or if the first bit is 1, returns:
       ('NEW', int_in_last_6_bits)
    """
    binary_str = f'{site_code:032b}'

    # If the first bit is 1 => "NEW"
    if int(binary_str[0]) == 1:
        return ('NEW', int(binary_str[-6:], 2))

    # Otherwise, skip the first 2 bits:
    #   bit[0] was 0, bit[1] we discard as well
    #   then parse region(17 bits), position(9 bits), strand(1 bit), damage(3 bits)
    #   total so far = 17 + 9 + 1 + 3 = 30 bits
    binary_str = binary_str[2:]  # skip 2 bits

    region = int(binary_str[:17], 2)
    pos = int(binary_str[17:26], 2)
    strand_bit = int(binary_str[26])
    strand = '+' if strand_bit == 0 else '-'
    dmg = int(binary_str[27:30], 2)

    return (region, pos, strand, dmg)


def pack_tuple_nochrom(tup):
    """
    Re-encode (region, pos, strand, damage) back into a 32-bit int
    matching the original "no-chrom" scheme.
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

    # Combine them exactly as decode_nochrom expects after skipping 2 bits
    # i.e. region + pos + strand + dmg (30 bits).
    # We'll store them in the lower 30 bits of the final 32-bit integer.
    # The top 2 bits are '00' for normal lines. '10...' for 'NEW' sentinel.
    bit_str = region_bits + pos_bits + strand_bit + dmg_bits
    int_val = int(bit_str, 2)
    return int_val


###############################################################################
#                         2) FASTA & INTERVALS (Optional)                     #
###############################################################################

def read_fasta(fasta_path: str):
    """
    Reads a FASTA file into a dict: {chrom: sequence}.
    Not strictly needed unless you are doing sequence checks.
    """
    genome = {}
    with open(fasta_path, 'r') as f:
        chrom = None
        seq = []
        for line in f:
            if line.startswith('>'):
                if chrom is not None:
                    genome[chrom] = ''.join(seq)
                chrom = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if chrom is not None:
            genome[chrom] = ''.join(seq)
    return genome


def hold_interval_indices(interval_file: str):
    """
    Parses a BED-like file with columns: chrom, start, ...
    Returns a list of (chrom, start), used for region indexing.
    """
    intervals = []
    with open(interval_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                intervals.append((chrom, start))
    return intervals


###############################################################################
#                       3) PARSING THE INPUT .BIN FILE                        #
###############################################################################

def parse_input_bin(file_path: str):
    """
    Single-pass read of the entire input .bin file, returning:
       data_dict[id] = np.ndarray of shape (N, 3) 
         storing (region, pos, strand_bit) for each line of that ID.

    We do NOT store damage here.
    We do NOT break on "NEW 0" -- we read until EOF, capturing all IDs.
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
            decoded = decode_nochrom(val)

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
                region_val, pos_val, strand_sym, _ = decoded
                strand_bit = 0 if strand_sym == '+' else 1
                block_rows.append((region_val, pos_val, strand_bit))

        # (End while) we already stored the last block above if any
    return data_dict


###############################################################################
#                          4) DAMAGE REDISTRIBUTION                           #
###############################################################################

def redistribute_uniform(array_3col: np.ndarray, total_dmg: int):
    """
    Distribute `total_dmg` across array_3col (N,3) => (region, pos, strand_bit)
    using a uniform multinomial.

    Returns an (M,4) array => (region, pos, strand_bit, damage) for M <= N
    (only sites with non-zero damage).
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
#          5) PERFORM SIMULATIONS (generating multiple run .bin files)        #
###############################################################################

def perform_simulations(
    file_path: str,
    index_path: str,
    runs: int,
    start_run_id: int,
    output_prefix: str
):
    """
    1) Parse the `.bin` file once => data_dict[id] = Nx3 array of (region,pos,strand_bit).
    2) Read the `index_path` (a table with columns for damage). We use dmg_total for each ID.
    3) For each run (1..runs):
       - Open `acc_run_full_<runID>.bin` for writing
       - For each ID in ascending order:
         = Write "NEW" sentinel => 2^31 + ID
         = Redistribute damage from index => Nx3 => Nx1 => Nx4
         = Write each row (region, pos, strand, damage) if damage>0
       - Write one final "NEW" sentinel => 2^31 + (max_id + 1)
    """

    if runs < 1:
        print("[WARNING] 'runs' must be >=1. Doing nothing.")
        return

    # 1) Parse input bin
    print("[INFO] Parsing input .bin file...")
    data_dict = parse_input_bin(file_path)
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
                # 'NEW' sentinel = top bit set + cur_id in the lower bits
                sentinel_code = (1 << 31) + cur_id
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
                    encoded = pack_tuple_nochrom((r_val, p_val, strand_char, dmg_val))
                    out_f.write(struct.pack('I', encoded))

            # Final sentinel => 'NEW' with (max_id + 1)
            final_sentinel = (1 << 31) + (max_id + 1)
            out_f.write(struct.pack('I', final_sentinel))

        print(f"   [INFO] Finished run {run_id} => {out_fname}")

    elapsed = time.time() - start_t
    print(f"[INFO] Completed {runs} run(s) in {elapsed:.2f} sec total.")


###############################################################################
#                           6) OPTIONAL DECODER CHECK                          #
###############################################################################

def decode_result(site_code: int):
    """
    Similar to decode_nochrom, but used on the final 'result' .bin structure,
    which typically sets bit[0] = 1 for 'NEW' or uses a different offset for damage bits.
    This is basically the same logic you had originally in decode_result.
    """
    bs = f"{site_code:032b}"
    if bs[0] == '1':
        return ("NEW", int(bs[-6:], 2))

    # skip the first bit
    bs = bs[1:]
    region = int(bs[:17], 2)
    pos    = int(bs[17:26], 2)
    strand_bit = int(bs[26])
    strand = '+' if strand_bit == 0 else '-'
    dmg = int(bs[27:31], 2)
    return (region, pos, strand, dmg)


def simple_decode_check(bin_path: str, index_path: str):
    """
    Minimal checker: reads the generated .bin file, sums the damage for each ID,
    and compares to the index file if the ID is < len(index).
    """
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
            dec = decode_result(val)
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
                dmg_val = dec[3]
                total_dmg += dmg_val


###############################################################################
#                           7) MAIN (CLI INTERFACE)                           #
###############################################################################

if __name__ == "__main__":
    """
    Example usage:
      python big_redistribute.py \
        --id 1 \
        --runs 10 \
        --input encodings/atac_4mers_sorted_nb.bin \
        --index ../../data/encodings/atac_nb_dmg_index.out \
        --output /usr/xtmp/bc301/sim_exp/

    After runs, you can optionally call `simple_decode_check` on one .bin file:
      python big_redistribute.py --decode /usr/xtmp/bc301/sim_exp/acc_run_full_1.bin \
        --index ../../data/encodings/atac_nb_dmg_index.out
    """
    parser = argparse.ArgumentParser(description="NumPy-based uniform redistribution of 4-mer damages.")
    parser.add_argument("--id", type=int, required=True,
                        help="Starting run ID offset (e.g., 1 => first run is '1', second run is '2', etc.)")
    parser.add_argument("--runs", type=int, default=10,
                        help="Number of runs to perform in one go.")
    parser.add_argument("--input", type=str, required=True,
                        help="Path to the input .bin file (4-mer data in no-chrom format).")
    parser.add_argument("--index", type=str, required=True,
                        help="Path to the .out or .csv with columns for damage (dmg_total).")
    parser.add_argument("--output", type=str, default=".",
                        help="Directory or prefix where output .bin files will go.")
    parser.add_argument("--decode", type=str, default="",
                        help="If provided, decode this file (instead of performing simulations).")

    args = parser.parse_args()

    if args.decode:
        # Just decode-check
        if not os.path.isfile(args.decode):
            print(f"[ERROR] decode file not found: {args.decode}")
        else:
            simple_decode_check(args.decode, args.index)
    else:
        # Perform simulations
        perform_simulations(
            file_path=args.input,
            index_path=args.index,
            runs=args.runs,
            start_run_id=args.id - 1,
            output_prefix=args.output
        )
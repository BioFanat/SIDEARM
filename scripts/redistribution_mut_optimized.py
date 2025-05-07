import argparse
import os
import struct
import time

import numpy as np
import pandas as pd

###############################################################################
#                           1) DECODING / ENCODING                            #
###############################################################################
def decode_mutation(site_code: int):
    """
    Decodes a 32-bit integer (site_code) into:
       ("NEW", new_id)  if the top bit is 1
    or
       (region, pos, damage)  if the top bit is 0

    Layout (for normal lines):
      - bit[0] = 0
      - region: next 16 bits
      - pos: next 9 bits
      - damage: next 6 bits
      ( total = 1 + 16 + 9 + 6 = 32 )
    """
    binary_str = f'{site_code:032b}'
    if binary_str[0] == '1':
        # "NEW" sentinel
        # lower 31 bits = ID
        new_id = int(binary_str[1:], 2)
        # print(binary_str, new_id)
        return ("NEW", new_id)
    else:
        # parse normal line
        # skip the first bit
        bits = binary_str[1:]  # now length 31
        region_val = int(bits[0:16], 2)
        pos_val    = int(bits[16:25], 2)
        dmg_val    = int(bits[25:31], 2)
        return (region_val, pos_val, dmg_val)

def encode_mutation(region: int, pos: int, damage: int):
    """
    Re-encodes (region, pos, damage) into a 32-bit integer with the top bit = 0.
    Layout:
      bit[0] = 0
      region (16 bits)
      pos (9 bits)
      damage (6 bits)
    """
    # Basic checks
    if not (0 <= region < 2**16):
        raise ValueError(f"[encode_mutation] region={region} out of 16-bit range.")
    if not (0 <= pos < 2**9):
        raise ValueError(f"[encode_mutation] pos={pos} out of 9-bit range.")
    if not (0 <= damage < 2**6):
        raise ValueError(f"[encode_mutation] damage={damage} out of 6-bit range.")

    region_bits = f"{region:016b}"
    pos_bits    = f"{pos:09b}"
    dmg_bits    = f"{damage:06b}"
    full_str    = "0" + region_bits + pos_bits + dmg_bits  # length = 32
    return int(full_str, 2)

def encode_new_id(new_id: int):
    """
    Encodes a 'NEW' sentinel with the highest bit = 1 
    and the lower 31 bits = new_id.
    """
    if not (0 <= new_id < 2**31):
        raise ValueError(f"[encode_new_id] new_id={new_id} out of 31-bit range.")
    return (1 << 31) | new_id


###############################################################################
#                       2) PARSE INPUT .BIN INTO MEMORY                       #
###############################################################################
def parse_input_bin(bin_path: str):
    """
    Reads the entire input .bin file (with "NEW" sentinel logic),
    returning a dict:
        data_dict[id] = np.array( [ [region, pos], [region, pos], ... ] )
    The ID is determined by the lower 31 bits of the sentinel line.
    We stop reading if we see a 'NEW' with ID=0 as an end sentinel (typical).
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
                
                break

            val = struct.unpack('I', chunk)[0]
            decoded = decode_mutation(val)

            if decoded[0] == "NEW":
                # store the previous block
                if current_id is not None and len(block_rows) > 0:
                    arr_np = np.array(block_rows, dtype=np.uint32)
                    data_dict[current_id] = arr_np

                new_id = decoded[1]
                if new_id == 16:
                    # typical stopping condition
                    break

                current_id = new_id
                block_rows = []
            else:
                # normal line => (region, pos, dmg)
                # but we ignore 'dmg' in the input for the re-distribution step
                region_val, pos_val, _ = decoded
                block_rows.append([region_val, pos_val])
        
        if current_id is not None and len(block_rows) > 0:
            arr_np = np.array(block_rows, dtype=np.uint32)
            data_dict[current_id] = arr_np

    return data_dict


###############################################################################
#                            3) READ INDEX FILE                               #
###############################################################################
def read_index_file(index_path: str):
    """
    Suppose the index file has one column (or more) 
    but the first/only relevant column is the total mutation count.

    Returns a 1D numpy array of length = number_of_rows
    such that index_array[id] = total_mutations.
    If your ID is the row number, this will match straightforwardly.
    """
    df = pd.read_csv(index_path, sep='\t', header=None)
    # If the first column is the total mutation count:
    # total_col = df[0].values
    # Or if the second or third column is the relevant one, adjust accordingly.
    total_col = df[0].values  # <— Adjust if necessary

    return total_col


###############################################################################
#                    4) UNIFORM REDISTRIBUTION OF MUTATIONS                   #
###############################################################################
def redistribute_uniform(array_2col: np.ndarray, total_mut: int):
    """
    Distributes 'total_mut' mutations uniformly across 'array_2col' (N,2),
    i.e. each row has probability 1/N of receiving each mutation.

    Returns an array of shape (M,3): [region, pos, damage], 
    only for rows that got a nonzero damage.
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

    out = np.column_stack([ array_2col[nz_mask], counts[nz_mask] ])
    return out


###############################################################################
#           5) PERFORM SIMULATIONS, WRITING MULTIPLE .BIN OUTPUTS            #
###############################################################################
def perform_simulations(
    input_bin: str,
    index_path: str,
    runs: int,
    start_run_id: int,
    output_prefix: str
):
    """
    - Parse the input .bin once -> data_dict[ID] = Nx2 array (region, pos)
    - Read the index_path -> array_of_mutations[ID] = total_mut
    - For i in [1..runs]:
       -> create an output .bin file (something like "mut_run_full_<runID>.bin")
       -> for each ID in ascending order:
            => write 'NEW' sentinel (1<<31) + ID
            => uniform redistribute total_mut from index array
            => for each site that got allocated >0, write (region,pos,damage)
       -> final sentinel (ID = max_id + 1)
    """

    if runs < 1:
        print("[WARNING] 'runs' must be >= 1. No simulations performed.")
        return

    # 1) parse the .bin
    print("[INFO] Parsing input bin...")
    data_dict = parse_input_bin(input_bin)
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
        # e.g. run_id = start_run_id * runs + (i + 1) 
        # or some other scheme you prefer
        run_id = (start_run_id-1) * runs + (i + 1)
        out_fname = os.path.join(output_prefix, f"mut_run_full_{run_id}.bin")

        with open(out_fname, 'wb') as out_f:
            for cur_id in all_ids:
                # write NEW sentinel
                sentinel_code = encode_new_id(cur_id)
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
                    region_val = row[0]
                    pos_val    = row[1]
                    dmg_val    = row[2]
                    encoded = encode_mutation(region_val, pos_val, dmg_val)
                    out_f.write(struct.pack('I', encoded))

            # final sentinel => NEW with (max_id+1) or 0 if you want to end
            final_sentinel = encode_new_id(max_id+1)
            out_f.write(struct.pack('I', final_sentinel))

        print(f"   [INFO] Finished run {i+1} => {out_fname}")

    elapsed = time.time() - start_t
    print(f"[INFO] Completed {runs} run(s) in {elapsed:.2f} sec total.")


###############################################################################
#                            6) OPTIONAL DECODER                              #
###############################################################################
def decode_result(bin_path: str, index_path: str):
    """
    A minimal check. We read the generated .bin file and sum the total damage 
    for each ID (only if ID < len(index_array)). Then compare to what's in the index.
    """
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
            dec = decode_mutation(val)
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
                # normal line => add damage
                dmg_val = dec[2]
                total_dmg_for_id += dmg_val


###############################################################################
#                               7) MAIN CLI                                   #
###############################################################################
if __name__ == "__main__":
    """
    Example usage:
      python mutation_redistribute.py \
        --input path/to/mutation_input.bin \
        --index path/to/mutation_index.tsv \
        --output /tmp/sim_output \
        --runs 10 \
        --start_id 1

    Then you get output files: /tmp/sim_output/mut_run_full_{1..10}.bin
    """
    parser = argparse.ArgumentParser(description="Uniform re-distribution of mutations (similar to the damage script).")
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
    parser.add_argument("--check", type=str, default="",
                        help="If provided, will run a decode check on that .bin instead of performing simulations.")

    args = parser.parse_args()

    if args.check:
        decode_result(args.check, args.index)
    else:
        perform_simulations(
            input_bin=args.input,
            index_path=args.index,
            runs=args.runs,
            start_run_id=args.start_id,
            output_prefix=args.output
        )
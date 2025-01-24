import struct
import time
import numpy as np
import pandas as pd
import argparse
import pickle
import gzip

def redistribute(lists, total_items):
    """
    Redistribute a total number of damage items across a list of 4-mers using pure random distribution.
    
    Args:
        lists: List of lists containing [region, position, strand, damage] information
        total_items: Total number of damages to redistribute
        
    Returns:
        List of lists with redistributed damage counts
    """
    # Extract columns from input lists
    region, pos, strand, dmg = [row[0] for row in lists], [row[1] for row in lists], \
                              [row[2] for row in lists], [float(row[3]) for row in lists]
    
    intersect = np.array([pos, dmg]).T
    num_destinations = len(intersect[:, 0])
    
    # Ensure total_items is an integer
    total_items = int(round(total_items))
    
    if num_destinations == 0:
        return []
    
    if num_destinations == 1:
        return [list(x) for x in zip(region, pos, strand, [total_items])]
    
    # Define uniform probabilities
    probabilities = np.ones(num_destinations) / num_destinations
    
    # Ensure probabilities sum to exactly 1
    probabilities = probabilities / np.sum(probabilities)
    
    # Randomly distribute all items at once
    destinations = np.random.choice(
        num_destinations,
        size=total_items,
        p=probabilities,
        replace=True
    )
    
    # Count items per destination using bincount
    counts = np.bincount(destinations, minlength=num_destinations)
    
    # Verification step
    assert np.sum(counts) == total_items, \
        f"Error: Distributed {np.sum(counts)} items instead of {total_items}"
    
    # Convert counts to list and combine with other data
    return [list(x) for x in zip(region, pos, strand, counts.tolist())]

def decode_nochrom(site_code):
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', int(string[-6:], 2))

    string = string[2:]

    region = int(string[:17], 2)
    pos = int(string[17:26], 2)

    strand = int(string[26])
    match strand:
        case 0:
            strand = '+'
        case 1:
            strand = '-'
    
    dmg = int(string[27:30], 2)

    return (region, pos, strand, dmg)

def decode_result(site_code):
    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', int(string[-6:], 2))

    string = string[1:]

    region = int(string[:17], 2)
    pos = int(string[17:26], 2)

    strand = int(string[26])
    match strand:
        case 0:
            strand = '+'
        case 1:
            strand = '-'
    
    dmg = int(string[27:31], 2)

    return (region, pos, strand, dmg)


def rev_complement(seq): # find the reverse complement for a sequence
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def pack_tuple_nochrom(tup):


    if not (0 <= tup[0] < 2**17):   # Check region fits in 16 bits
        raise ValueError(f"Region {tup[0]} too large for 16 bits")
    if not (0 <= tup[1] < 2**9):   # Check position fits in 9 bits
        raise ValueError(f"Position {tup[1]} too large for 9 bits")
    if not (0 <= tup[3] < 2**4):   # Check damage fits in 4 bits
        raise ValueError(f"Damage {tup[3]} too large for 4 bits")

    region = f'{tup[0]:017b}'

    pos = f'{tup[1]:09b}'

    strand = ''

    match tup[2]:
        case '+':
            strand = '0'
        case '-':
            strand = '1'

    dmg = f'{tup[3]:04b}'

    return int(region + pos + strand + dmg, 2)

def hold_interval_indices(interval_file):
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            
            start = int(cur[1])
            intervals.append((cur[0], start))
    
    return intervals


def convert_id(id, bpde=False):
    # convert id to base 4
    id_str = np.base_repr(id, base=4).zfill(3 - int(bpde))

    ret = []
    match id_str[0]:
        case '0':
            ret.append('A')
        case '1':
            ret.append('C')
        case '2':
            ret.append('G')
        case '3':
            ret.append('T')


    if not bpde:
        match id_str[1]:
            case '0':
                ret.append('CC')
            case '1':
                ret.append('CT')
            case '2':
                ret.append('TC')
            case '3':
                ret.append('TT')

        match id_str[2]:
            case '0':
                ret.append('A')
            case '1':
                ret.append('C')
            case '2':
                ret.append('G')
            case '3':
                ret.append('T')
    
    else:
        ret.append('G')

        match id_str[1]:
            case '0':
                ret.append('A')
            case '1':
                ret.append('C')
            case '2':
                ret.append('G')
            case '3':
                ret.append('T')
    return ''.join(ret)

def perform_run(file_path, index, id, genome, intervals): #bin_file is a file stream
    bin_file = open(file_path, 'rb')

    run_file = open(f'/usr/xtmp/bc301/sim_runtime/acc_run_full_{id}.bin', 'wb')
    # unpacked_run = open(f'/usr/xtmp/bc301/sim_exp/acc_run_full_{id}.bed', 'w')
    # pkl_file = open(f'/usr/xtmp/bc301/sim_exp/acc_run_full_{id}.pkl', 'wb')
    # # run_file = open(f'encodings/simulated_atac_run_{id}.bin', 'wb')

    data = bin_file.read(4)
    source = decode_nochrom(struct.unpack('I', data)[0])
    
    current_id = source[1]
    arr = []

    run_file.write(struct.pack('I', 2**31 + current_id))

    while data:
        data = bin_file.read(4)

        source = decode_nochrom(struct.unpack('I', data)[0])
        
        if source[0] == "NEW":

            if len(arr) > 0:
                arr = redistribute(arr, index[current_id])
                for a in arr:
                    # chrom, start = intervals[a[0]]
                    # seq = genome[chrom][start + a[1] - 1:start + a[1] + 2]

                    # strand = a[2]
                    # if strand == '+' and seq.upper() != convert_id(current_id, bpde=True) or strand == '-' and seq.upper() != rev_complement(convert_id(current_id, bpde=True)):
                    #     print(f"UNMATCHED for {current_id} - SHOULD BE {convert_id(current_id, bpde=True)}, GOT {seq}")
                    
                    if a[3] > 0:
                        run_file.write(struct.pack('I', pack_tuple_nochrom(a)))
                        # unpacked_run.write(f"{chrom}\t{start + a[1]}\t{start + a[1] + 2}\t{a[2]}\t{a[3]}\n")
                    # Also write to pickled format for intermediate compression
                    # pkl_file.write(pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL))
                    # with gzip.open(f'/usr/xtmp/bc301/sim_exp/acc_run_full_{id}.pkl.gz', 'ab') as f:
                    #     pickle.dump(a, f)
            if source[1] == 0:
                break

            arr = []
            current_id = source[1]
            
            run_file.write(struct.pack('I', 2**31 + current_id))
            
        else:
            arr.append(source)
        
    run_file.write(struct.pack('I', 2**31 + current_id + 1))
    run_file.close()
    bin_file.close()
    # unpacked_run.close()
    # pkl_file.close()

def read_fasta(fasta_file): # reads genomic sequences - in .fasta format
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

def perform_simulations(file_path, runs, start):

    hg19 = read_fasta("../../data/genome/hg19.fa")

    intervals = hold_interval_indices("../../data/raw/atac_150bp_intervals_merged.bed")

    acc_4mer = pd.read_csv('../../data/encodings/atac_nb_dmg_index.out', sep="\t", header=None, names=["dmg_plus", "dmg_minus", "counts_plus", "counts_minus"])
    acc_4mer['dmg_total'] = acc_4mer['dmg_plus'] + acc_4mer['dmg_minus']
    index = acc_4mer['dmg_total'].values
    
    start_time = time.time()

    for i in range(1, runs + 1):
        perform_run(file_path, index, start*runs + i, hg19, intervals)
    
    print(f"whole operation took {time.time() - start_time:.2f} seconds")

def simple_decode(file_path, genome, intervals):
    bin_file = open(file_path, 'rb')
    total_dmg = 0
    max_dmg = 0    

    acc_4mer = pd.read_csv('../../data/encodings/atac_nb_dmg_index.out', sep="\t", header=None, names=["dmg_plus", "dmg_minus", "counts_plus", "counts_minus"])
    acc_4mer['dmg_total'] = acc_4mer['dmg_plus'] + acc_4mer['dmg_minus']
    index = acc_4mer['dmg_total'].values

    current_id = 0

    total_dmg, max_dmg = 0, 0

    while True:
        data = bin_file.read(4)
        if len(data) < 4:  # Check if we got a full 4 bytes
            break
            
        source = decode_result(struct.unpack('I', data)[0])
        

        if source[0] == "NEW":
            # if current_id > 0 and current_id <= len(index) and total_dmg != index[current_id-1]:
            #     # print(current_id-1, total_dmg, max_dmg, index[current_id-1])
            #     print(f"STARTING NEW ID {convert_id(current_id)} or {rev_complement(convert_id(current_id))}")
            print(total_dmg, max_dmg, current_id-1)
    
            total_dmg, max_dmg = 0, 0
            current_id += 1

        else:    
            chrom, start = intervals[source[0]]
            seq = genome[chrom][start + source[1] - 1:start + source[1] + 2] 
            # print("we're using", seq, "and a strand of", source[2])
            total_dmg += source[-1]
            max_dmg = max(max_dmg, source[-1])
    
    bin_file.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A script that accepts command-line inputs")

    parser.add_argument("id", type=int, help="Starting Run Id")
    args = parser.parse_args()

    # perform_simulations('encodings/acc_4mers_atac_seq_sorted.bin', 100, args.id - 1)

    perform_simulations('../../data/encodings/atac_4mers_sorted_nb.bin', 1000, args.id - 1)

    # hg38 = read_fasta("../../data/genome/hg38.fa")
    # intervals = hold_interval_indices("../../data/raw/a549_regions_merged.bed")

    # simple_decode('/usr/xtmp/bc301/sim_bpde_nb/acc_run_1.bin', hg38, intervals)

    # hg19 = read_fasta("../../data/genome/hg19.fa")
    # intervals = hold_interval_indices("../../data/raw/atac_150bp_intervals_merged.bed")

    # simple_decode('/usr/xtmp/bc301/sim_uv_cpd_full/acc_run_full_3001.bin', hg19, intervals)

    # print(f"Decoding run {args.id}")
    # simple_decode(f"encodings/simulated_atac_run_{args.id}.bin")
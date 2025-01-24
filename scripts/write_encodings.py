from genome_utils.general_utils import read_fasta, rev_complement
import pandas as pd
import struct
import numpy as np
import time
import argparse

def match_C_3mer(seq):
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

def sliding_window(key, seq, categories, region_id, mode):

    match mode:
        case 'C_mut':
            length = 3
            match_func = match_C_3mer
        case 'CPD_mut':
            length = 4
            match_func = match_CPD_4mer

    if len(seq) < length:
        return categories
    
    split = key.split(":")
    chrom = split[-2]
    start = int(split[-1].split("-")[0])
    
    for i in range(len(seq) - length + 1):
        cur = seq[i:i+length].upper()
        index = match_func(cur)
        if index < 16:
            result_tuple = (region_id, i + 1, 0) #(chrom, i + 1, site_strand, region_id, 0)
            categories[index].append(result_tuple)
    
    return categories

def compress_3mer(tup, region_bits=16, pos_bits=9, strand_bits=0):
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

def write_compression(index, output='../data/encodings/uv_fibroblasts_3mer_packed.bin'):
    bin_file = open(output, 'wb')
    index_id = 0
    # max_region = 0
    for seq_id in index:
        bin_file.write(struct.pack('I', 2**31 + index_id))

        for site in seq_id:
            
            bin_file.write(struct.pack('I', compress_3mer(site)))
            # max_region = max(max_region, decompress_3mer(compress_3mer(site))[0])
           
        index_id += 1
    
    bin_file.write(struct.pack('I', 2**31 + index_id))

def decompress_3mer(site_code, region_bits=16, pos_bits=9, strand_bits=0):
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


def generate_categories(genome, mutations, output):
    counts = [0] * 16
    with open(mutations) as f:
        for line in f:
            mut = line.strip().split("\t")
            if mut[0] in genome:
                seq = genome[mut[0]][int(mut[1])-1:int(mut[1])+2]
                val = match_C_3mer(seq)
                if val < 16:
                    counts[val] += 1
                # if seq[1].upper() != mut[-3]:
                #     print(f'{mut[1]}\t{seq}\t{mut[-3]}') 

    out = open(output, 'w')
    i = 0
    for ind in counts:
        # seq = "\t".join(decode_3mer(i))
        out.write(f'{ind}\n')      
        i += 1
    
    out.close()

def perform_run(file_path, index, id): #uses MUTATIONS instead
    bin_file = open(file_path, 'rb')

    run_file = open(f'/usr/xtmp/bc301/sim_data_skin_mut_atac/acc_run_{id}.bin', 'wb')
    #run_file = open(f'encodings/simulated_mut_{id}.bin', 'wb')

    data = bin_file.read(4)
    source = decompress_3mer(struct.unpack('I', data)[0])
    
    current_id = source[1]

    #print(current_id, "-", index[current_id])
    arr = []

    run_file.write(struct.pack('I', 2**31 + current_id))
    #dmg_sum = 0
    max_dmg = 0
    max_region = 0

    while data:
        
        data = bin_file.read(4)

        source = decompress_3mer(struct.unpack('I', data)[0])
        
        if source[0] == "NEW":
            print(current_id, "-", index[current_id])
            
            dmg_sum = 0
            if len(arr) > 0:
                arr = redistribute(arr, index[current_id])
                for a in arr:
                    run_file.write(struct.pack('I', compress_3mer(a)))
                    dmg_sum += a[-1]
                    max_dmg = max(max_dmg, a[-1])
                    max_region = max(max_region, a[0])
                    
    
            if source[1] == 0:
                break
            
            if dmg_sum != index[current_id]:
                print(current_id, "-", index[current_id])

            arr = []
            current_id = source[1]
            
            run_file.write(struct.pack('I', 2**31 + current_id))
        else:
            arr.append(source)
        
    run_file.write(struct.pack('I', 2**31 + current_id + 1))
    run_file.close()
    bin_file.close()
    
    print(max_dmg)
    print(max_region)

def redistribute(lists, total_items): # redistribute the total damage within a sequence to its components

    region, pos, dmg = [row[0] for row in lists], [row[1] for row in lists], [float(row[2]) for row in lists]
    
    intersect = np.array([pos, dmg]).T

    num_destinations = len(intersect[:, 0])

    # Define the probabilities of selecting each destination
    probabilities = np.ones(num_destinations) / num_destinations

    # Distribute the items randomly

    destinations = np.random.choice(num_destinations, int(total_items), p=probabilities)

    # Count the number of items per destination
    counts = np.bincount(destinations, minlength=num_destinations)

    # make counts all floats
    counts = counts.astype(int).tolist()

    # Print the final distribution
    return [list(x) for x in zip(region, pos, counts)]

def perform_simulations(file_path, runs, start):
    # FOR MUTATIONS sums = pd.read_csv('mutations_by_3mer_aggregated_acc_cor.out', header=None)[0].values
    sums = pd.read_csv('../../data/mutation_skin_3mer_accessible.out', header=None)[0].values # FOR MUTAGENIC
    
    #perform_run(file_path, sums, 1)
    
    start_time = time.time()

    for i in range(1, runs + 1):
        perform_run(file_path, sums, start*runs + i)
    
    print(f"whole operation took {time.time() - start_time:.2f} seconds")

def simple_decode(file_path):
    bin_file = open(file_path, 'rb')
    total_dmg = 0
    max_dmg = 0    

    while True:
        data = bin_file.read(4)
        if len(data) < 4:  # Check if we got a full 4 bytes
            break
            
        source = decompress_3mer(struct.unpack('I', data)[0])

        if source[0] == "NEW":
            print(total_dmg, max_dmg)
            total_dmg, max_dmg = 0, 0
        else:    
            total_dmg += source[-1]
            max_dmg = max(max_dmg, source[-1])
    
    print(total_dmg, max_dmg)
    bin_file.close()


if __name__ == "__main__":

    # genome = read_fasta("../../data/genome/hg19.fa")
    # generate_categories(genome=genome, mutations='../../data/raw/atac_mutations_transitions_C_only.bed', output='../../data/mutation_skin_3mer_accessible.out')

    # regions = read_fasta("../../data/raw/atac_150bp_intervals_merged.fa.out")

    # pos = [[] for i in range(16)]
    # region_id = 0

    # for key in list(regions.keys()):
    #     pos = sliding_window(key, regions[key], pos, region_id, 'C_mut')
    #     region_id += 1
    
    # print(region_id)

    # write_compression(pos, '../../data/encodings/skin_mutations_packed_atac.bin')


    # print("compression successful")

    parser = argparse.ArgumentParser(description="A script that accepts command-line inputs")

    parser.add_argument("id", type=int, help="Starting Run Id")
    args = parser.parse_args()

    perform_simulations('../../data/encodings/skin_mutations_packed_atac.bin', 100, args.id-1)

    # simple_decode('/usr/xtmp/bc301/sim_data_skin_mut_atac/acc_run_1.bin')

    print("simulation successful")

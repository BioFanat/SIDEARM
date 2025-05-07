import pickle
import struct
import numpy as np

def collapse_dmg_sites(input_file): # reads CPDSeq 2.0 damages - in .bed format
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

    #f = open(output_file, 'w')

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

def write_compressed_damages(sites, output):
    file = open(output, 'w')
    for site in sites:
        file.write("\t".join([str(word) for word in site]) + "\n")
    
    file.close()

def convert_id(id):
    # convert id to base 4
    id_str = np.base_repr(id, base=4).zfill(3)

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
    return ''.join(ret)

def rev_complement(seq): # find the reverse complement for a sequence
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

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

def sliding_window(key, seq, categories, region_id, intervals, genome):

    length = 4
    if len(seq) < length:
        return categories
    
    split = key.split(":")
    chrom = split[0]
    start = int(split[-1].split("-")[0])

    # test_chrom, test_start = intervals[region_id]
    # if chrom != test_chrom or start != test_start:
    #     print(f"MISMATCH for {region_id}: {chrom}:{start} != {test_chrom}:{test_start}")


    
    for i in range(len(seq) - length + 1):
        cur = seq[i:i+length].upper()
        index, site_strand = match_cpd_4mer(cur)
        if index < 64:
    
            core = cur[1]
            if site_strand == "+" and core != "G":
                # This is physically impossible - can't have CPD at G/A on + strand
                print(f"Warning: Invalid BPDE site found: {cur} on + strand") #- matched to {convert_id(index)}, {site_strand}")
            if site_strand == "-" and core != "C":
                # This is physically impossible - can't have CPD at C/T on - strand
                print(f"Warning: Invalid BPDE site found: {cur} on - strand") #- matched to {convert_id(index)}, {site_strand}")

            # print(f"4mer: {cur}")
            # print(f"Index: {index}")
            # print(f"Strand: {site_strand}")

            #result_tuple = (chrom, i + start + 1, site_strand, 0)
            result_tuple = (region_id, i + 1, site_strand, 0) #(chrom, i + 1, site_strand, region_id, 0)


            encoded = pack_tuple_nochrom(result_tuple)  # Note: should use pack_tuple_nochrom since tuple format changed
            test_seq, test_chrom, test_pos = decode_to_seq(encoded, intervals, genome)

            #Now comparing the actual genomic coordinates
            # if (f'{chrom}' != test_chrom) or (start + i + 1 != test_pos):
            #     print(f"MISMATCH for {region_id}: {chrom}:{start + i + 1} != {test_chrom}:{test_pos}")
            
            categories[index].append(result_tuple)
    
    return categories

def map_collapsed_cpd(chrom, seq, dmg, pos):
    index, site_strand = match_cpd_4mer(seq)
    if index < 64:
        match site_strand:
            case "+":
                pos[index][0] += dmg
            case "-":
                pos[index][1] += dmg
        
    return pos

def map_collapsed_bpde(chrom, seq, dmg, pos):
    index, site_strand = match_bpde_3mer(seq)
    if index < 16:
        match site_strand:
            case "+":
                pos[index][0] += dmg
            case "-":
                pos[index][1] += dmg
        
    return pos

def is_valid_chrom(chrom):
    diff = chrom[3:]
    return chrom[:3] == 'chr' and diff in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

def pack_tuple_nochrom(tup):
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

def decode_to_seq(site_code, interval_map, genome):

    string = f'{site_code:032b}'

    if int(string[0]) == 1:
        return ('NEW', f'{int(string[-6:-4], 2)}{int(string[-4:-2], 2)}{int(string[-2:], 2)}')

    string = string[2:]

    region = int(string[:17], 2)
    pos = int(string[17:26], 2)

    chrom, start = interval_map[region]
    pos += start

    test_seq = genome[chrom][pos-1:pos+3]

    # strand = int(string[26])
    # match strand:
    #     case 0:
    #         strand = '+'
    #     case 1:
    #         strand = '-'
    
    # dmg = int(string[27:30], 2)

    return test_seq, chrom, pos

def hold_interval_indices(interval_file):
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            
            start = int(cur[1])
            intervals.append((cur[0], start))
    
    return intervals


if __name__ == "__main__":

    # sites_plus = collapse_dmg_sites("/usr/xtmp/hiw4/CPDdata/XPC_12J_NakedDNA_S22_CPD_1bp_sorted_plus_reordered.bed")
    # sites_minus = collapse_dmg_sites("/usr/xtmp/hiw4/CPDdata/XPC_12J_NakedDNA_S22_CPD_1bp_sorted_minus_reordered.bed")
    
    # write_compressed_damages(sites_plus, "collapsed_naked_plus.bed")
    # write_compressed_damages(sites_minus, "collapsed_naked_minus.bed")
    hg19 = read_fasta("../../data/genome/hg19.fa")

    # intervals = hold_interval_indices("../../data/raw/atac_150bp_intervals_merged.bed")

    # regions = read_fasta("../../data/raw/atac_150bp_regions_merged.fa.out") #read_fasta("acc_regions/acc_intervals.fa.out") #read_fasta("acc_regions/acc_intervals_consolidated.fa.out")
    # # # # #print(list(regions.keys()))

    # pos = [[] for i in range(64)]
    # region_id = 0
    # max_region_id = 0

    # # #curr_chrom = ""
    # max_len = 0
    # for key in list(regions.keys()):
    #     #new_chrom = key.split(":")[0]

    #     # if new_chrom != curr_chrom:
    #     #     curr_chrom = new_chrom
    #     #     region_id = 0

    #     max_region_id = max(region_id, max_region_id)
    #     max_len = max(max_len, len(regions[key]))
    #     pos = sliding_window(key, regions[key], pos, region_id, intervals, hg38)
    #     region_id += 1
    
    # print(max_region_id)
    # print(max_len)

    # bin_file = open("../../data/encodings/atac_4mers_ultracompact.bin", "wb")
    # index_id = 0
    # bytes_written = 0  # Track total bytes written
    # site_by_index = [[0,0] for i in range(64)]
        
    # for index in pos:
    #     # Write marker
    #     bin_file.write(struct.pack('I', 2**31 + index_id))
    #     bytes_written += 4
        
        
    #     # Write sites
    #     for site in index:
    #         site_chrom, site_position = site[0], site[1]
    #         if site[3] > 0:

    #             packed = pack_tuple_nochrom(site)  
    #             test_seq, test_chrom, test_pos = decode_to_seq(packed, intervals, hg38)


    #             # if site_chrom != test_chrom or site_position != test_pos:
    #             #     print(f"MISMATCH: {site_chrom} != {test_chrom} or {site_position} != {test_pos}")

    #             #add sites based on strand

    #             if site[2] == '+':
    #                 site_by_index[index_id][0] += 1
    #                 if test_seq.upper() != convert_id(index_id):
    #                     print(f'Invalid CPD match: {test_seq} != {convert_id(index_id)}')
    #             else:
    #                 site_by_index[index_id][1] += 1
    #                 if test_seq.upper() != rev_complement(convert_id(index_id)):
    #                     print(f'Invalid CPD match: {test_seq} != {rev_complement(convert_id(index_id))}')

    #             bin_file.write(struct.pack('I', packed))
    #             bytes_written += 4
        
    #     index_id += 1
    
    # # Write final marker
    # bin_file.write(struct.pack('I', 2**31 + index_id))
    # bytes_written += 4

    # bin_file.close()

    # print("compression successful")

    dmg_index = [[0,0] for i in range(64)]

    with open('../../data/damages/atac_naked_plus.bed') as f:
        for line in f:
            arr = line.split("\t")
            chrom = arr[0]
            start = int(arr[1])
            end = int(arr[2])
            dmg = int(arr[3])

            if is_valid_chrom(chrom):
                seq = hg19[chrom][start-1:end+2]

                if len(seq) != 4:
                    print(line)

                dmg_index = map_collapsed_cpd(chrom, seq, dmg, dmg_index)
    
    with open('../../data/damages/atac_naked_minus.bed') as f:
        for line in f:
            arr = line.split("\t")
            chrom = arr[0]
            start = int(arr[1])
            end = int(arr[2])
            dmg = int(arr[3])
            
            if is_valid_chrom(chrom):
                seq = hg19[chrom][start-1:end+2]

                if len(seq) != 4:
                    print(line)

                dmg_index = map_collapsed_cpd(chrom, seq, dmg, dmg_index)
    
    file = open("../../data/encodings/atac_naked_dmg_index.out", "w")

    index_id = 0

    for ind in dmg_index:
        ind = [str(el) for el in ind]
        file.write("\t".join(ind) + "\n")
        index_id += 1
    
    file.close()
    print("compression successful")




    
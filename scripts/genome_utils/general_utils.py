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

def collapse_mut_sites(mut_path, output): # reads CPDSeq 2.0 damages - in .bed format
    

    out = open(output, 'w')
    #closed = True

    with open(mut_path) as infile:
        header = infile.readline()
        
        current_chrom, current_pos, current_orig, current_new = infile.readline().strip().split(",")
        current_pos = int(current_pos)
        current_chrom = f'chr{current_chrom}'
        count = 1

        for line in infile:
            chrom, pos, orig, new = line.strip().split(",")
            chrom = f'chr{chrom}'
            pos = int(pos)
            print(chrom, pos, current_chrom, current_pos, chrom == current_chrom, pos == current_pos, count)

            if chrom == current_chrom and pos == current_pos:
                count += 1

            else:
                out.write(f'{current_chrom}\t{current_pos-1}\t{current_pos-1}\t.\t{current_orig}\t{current_new}\t{count}\n')
                count = 1
            
            current_chrom, current_pos, current_orig, current_new = chrom, pos, orig, new

        out.write(f'{current_chrom}\t{current_pos-1}\t{current_pos}\t.\t{current_orig}\t{current_new}\t{count}\n')
                                        
    out.close()

def rev_complement(seq): # find the reverse complement for a sequence
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def is_valid_chrom(chrom):
    diff = chrom[3:]
    return chrom[:3] == 'chr' and diff in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']


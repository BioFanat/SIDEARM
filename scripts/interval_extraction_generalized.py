#!/usr/bin/env python3
"""
Generalized interval extraction for SIDEARM.

This script extracts k-mers from genomic intervals and encodes them
for efficient storage and processing. It supports different damage types
and k-mer configurations.
"""

import argparse
import pickle
import struct
import os
import numpy as np
from typing import Dict, List, Tuple, Any, Optional

# Import our new configuration system
from kmer_config import KmerConfig, get_config_by_name, get_default_config

def read_fasta(fasta_file: str) -> Dict[str, str]:
    """
    Read genomic sequences in FASTA format.
    
    Args:
        fasta_file: Path to the FASTA file
        
    Returns:
        Dictionary mapping chromosome names to sequences
    """
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
        if chrom is not None:
            genome[chrom] = seq
    return genome

def collapse_dmg_sites(input_file: str) -> List[List[Any]]:
    """
    Read CPDSeq 2.0 damages from a BED format file.
    
    Args:
        input_file: Path to the input BED file
        
    Returns:
        List of consolidated damage sites
    """
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

def write_compressed_damages(sites: List[List[Any]], output: str) -> None:
    """
    Write compressed damages to a file.
    
    Args:
        sites: List of damage sites
        output: Path to the output file
    """
    with open(output, 'w') as file:
        for site in sites:
            file.write("\t".join([str(word) for word in site]) + "\n")

def hold_interval_indices(interval_file: str) -> List[Tuple[str, int]]:
    """
    Read interval indices from a file.
    
    Args:
        interval_file: Path to the interval file
        
    Returns:
        List of tuples (chromosome, start_position)
    """
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            
            start = int(cur[1])
            intervals.append((cur[0], start))
    
    return intervals

def is_valid_chrom(chrom: str) -> bool:
    """
    Check if a chromosome name is valid.
    
    Args:
        chrom: Chromosome name
        
    Returns:
        True if valid, False otherwise
    """
    diff = chrom[3:]
    return chrom[:3] == 'chr' and diff in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

def sliding_window(key: str, seq: str, categories: List[List[Tuple]], region_id: int, 
                   intervals: List[Tuple[str, int]], genome: Dict[str, str], 
                   config: KmerConfig) -> List[List[Tuple]]:
    """
    Apply a sliding window to extract k-mers from a sequence.
    
    Args:
        key: Sequence key (e.g., "chr1:1000-2000")
        seq: DNA sequence
        categories: List of categories for k-mers
        region_id: Region ID
        intervals: List of interval indices
        genome: Genome dictionary
        config: KmerConfig object for k-mer matching
        
    Returns:
        Updated categories list
    """
    kmer_length = config.kmer_length
    if len(seq) < kmer_length:
        return categories
    
    split = key.split(":")
    chrom = split[0]
    start = int(split[-1].split("-")[0])
    
    for i in range(len(seq) - kmer_length + 1):
        cur = seq[i:i+kmer_length].upper()
        index, site_strand = config.match_kmer(cur)
        
        # Skip invalid k-mers
        if index >= config._get_invalid_kmer_code():
            continue
        
        # Create site tuple with region_id, position, strand, damage
        result_tuple = (region_id, i + 1, site_strand, 0)
        
        # Use the config to encode/decode for validation
        encoded = config.encode_site(region_id, i + 1, site_strand, 0)
        decoded = config.decode_site(encoded)
        
        # Skip if encoding/decoding didn't work correctly
        if decoded[0] != 'SITE' or decoded[1] != region_id or decoded[2] != i + 1:
            print(f"Warning: Encoding/decoding mismatch for {result_tuple} -> {decoded}")
            continue
            
        # Add to categories
        categories[index].append(result_tuple)
    
    return categories

def map_collapsed_dmg(chrom: str, seq: str, dmg: int, pos: List[List[int]], 
                       config: KmerConfig) -> List[List[int]]:
    """
    Map collapsed damage sites.
    
    Args:
        chrom: Chromosome name
        seq: DNA sequence
        dmg: Damage value
        pos: Position list
        config: KmerConfig object for k-mer matching
        
    Returns:
        Updated position list
    """
    index, site_strand = config.match_kmer(seq)
    
    # Skip invalid k-mers
    if index >= config._get_invalid_kmer_code():
        return pos
    
    # Update damage counts based on strand
    if site_strand == "+":
        pos[index][0] += dmg
    elif site_strand == "-":
        pos[index][1] += dmg
    
    return pos

def map_collapsed_dmg_generic(chrom: str, seq: str, dmg: int, pos: List[List[int]], 
                              config: KmerConfig) -> List[List[int]]:
    """
    Generic function to map collapsed damage sites for any k-mer type.
    
    Args:
        chrom: Chromosome name
        seq: DNA sequence
        dmg: Damage value
        pos: Position list
        config: KmerConfig object for k-mer matching
        
    Returns:
        Updated position list
    """
    index, site_strand = config.match_kmer(seq)
    
    # Check for valid k-mer
    if index >= config._get_invalid_kmer_code():
        return pos
    
    # Ensure index is within range
    if index >= len(pos):
        # Dynamically expand pos list if needed
        pos.extend([[0, 0] for _ in range(index - len(pos) + 1)])
    
    # Update damage counts based on strand
    if site_strand == "+":
        pos[index][0] += dmg
    elif site_strand == "-":
        pos[index][1] += dmg
    
    return pos

def process_input_files(fasta_path: str, plus_path: str, minus_path: str, intervals_path: str,
                        output_bin: str, output_index: str, damage_type: str = 'CPD') -> None:
    """
    Process input files to extract and encode k-mers.
    
    Args:
        fasta_path: Path to the reference genome FASTA file
        plus_path: Path to the plus strand damage file
        minus_path: Path to the minus strand damage file
        intervals_path: Path to the intervals file
        output_bin: Path to the output binary file
        output_index: Path to the output index file
        damage_type: Damage type (default: 'CPD')
    """
    # Get the configuration for this damage type
    try:
        config = get_config_by_name(damage_type)
    except ValueError:
        print(f"Warning: Unknown damage type '{damage_type}'. Using default (CPD).")
        config = get_default_config()
    
    print(f"Using configuration for {config.name} with kmer length {config.kmer_length}")
    
    # Load genome
    print("Loading genome...")
    genome = read_fasta(fasta_path)
    
    # Load intervals
    print("Loading intervals...")
    intervals = hold_interval_indices(intervals_path)
    
    # Determine max categories based on k-mer length
    max_categories = config._get_invalid_kmer_code()
    print(f"Using {max_categories} categories for {config.name}")
    
    # Initialize damage index
    dmg_index = [[0, 0] for _ in range(max_categories)]
    
    # Process plus strand damages
    print(f"Processing plus strand damages from {plus_path}...")
    with open(plus_path) as f:
        for line in f:
            arr = line.strip().split("\t")
            if len(arr) < 4:
                continue
                
            chrom = arr[0]
            start = int(arr[1])
            end = int(arr[2])
            dmg = int(arr[3])
            
            if is_valid_chrom(chrom):
                # Extract sequence with appropriate context for k-mer
                context_size = config.kmer_length - 1
                seq = genome[chrom][start-1:end+context_size]
                
                if len(seq) != config.kmer_length:
                    print(f"Warning: Invalid sequence length {len(seq)} at {chrom}:{start}-{end}")
                    continue
                
                dmg_index = map_collapsed_dmg_generic(chrom, seq, dmg, dmg_index, config)
    
    # Process minus strand damages
    print(f"Processing minus strand damages from {minus_path}...")
    with open(minus_path) as f:
        for line in f:
            arr = line.strip().split("\t")
            if len(arr) < 4:
                continue
                
            chrom = arr[0]
            start = int(arr[1])
            end = int(arr[2])
            dmg = int(arr[3])
            
            if is_valid_chrom(chrom):
                # Extract sequence with appropriate context for k-mer
                context_size = config.kmer_length - 1
                seq = genome[chrom][start-1:end+context_size]
                
                if len(seq) != config.kmer_length:
                    print(f"Warning: Invalid sequence length {len(seq)} at {chrom}:{start}-{end}")
                    continue
                
                dmg_index = map_collapsed_dmg_generic(chrom, seq, dmg, dmg_index, config)
    
    # Write damage index
    print(f"Writing damage index to {output_index}...")
    with open(output_index, "w") as file:
        for idx, ind in enumerate(dmg_index):
            ind = [str(el) for el in ind]
            file.write("\t".join(ind) + "\n")
    
    print("Damage index written successfully")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Generalized interval extraction for SIDEARM")
    parser.add_argument("--genome", required=True, help="Path to the reference genome FASTA file")
    parser.add_argument("--plus", required=True, help="Path to the plus strand damage file")
    parser.add_argument("--minus", required=True, help="Path to the minus strand damage file")
    parser.add_argument("--intervals", required=True, help="Path to the intervals file")
    parser.add_argument("--output-bin", required=True, help="Path to the output binary file")
    parser.add_argument("--output-index", required=True, help="Path to the output index file")
    parser.add_argument("--damage-type", default="CPD", help="Damage type (CPD, BPDE, or custom)")
    
    args = parser.parse_args()
    
    process_input_files(
        fasta_path=args.genome,
        plus_path=args.plus,
        minus_path=args.minus,
        intervals_path=args.intervals,
        output_bin=args.output_bin,
        output_index=args.output_index,
        damage_type=args.damage_type
    )

if __name__ == "__main__":
    main()
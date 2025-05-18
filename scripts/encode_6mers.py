#!/usr/bin/env python3
"""
Encode 6-mers from genomic intervals.

This script demonstrates how to extract and encode 6-mers from genomic sequences
using the generalized k-mer configuration system.
"""

import os
import struct
import argparse
from typing import Dict, List, Tuple, Any

import numpy as np

from kmer_config import MER6_CONFIG, get_config_by_name

def read_fasta(fasta_file: str) -> Dict[str, str]:
    """Read genomic sequences from a FASTA file."""
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

def hold_interval_indices(interval_file: str) -> List[Tuple[str, int]]:
    """Read interval indices from a file."""
    intervals = []
    with open(interval_file) as f:
        for line in f:
            cur = line.strip().split("\t")
            start = int(cur[1])
            intervals.append((cur[0], start))
    return intervals

def extract_6mers(regions: Dict[str, str], intervals: List[Tuple[str, int]], 
                  genome: Dict[str, str], output_bin: str, output_index: str):
    """
    Extract and encode 6-mers from genomic regions.
    
    Args:
        regions: Dictionary mapping region keys to sequences
        intervals: List of (chromosome, start) tuples
        genome: Dictionary mapping chromosome names to sequences
        output_bin: Path to output binary file
        output_index: Path to output index file
    """
    # Get the 6-mer configuration
    config = MER6_CONFIG
    
    # Maximum 6-mer code (4^6 = 4096 possible 6-mers)
    max_codes = 4096
    
    # Initialize categories
    categories = [[] for _ in range(max_codes)]
    
    # Process regions
    region_id = 0
    for key in regions.keys():
        seq = regions[key]
        
        if len(seq) < config.kmer_length:
            print(f"Warning: Region {key} is too short for 6-mers")
            region_id += 1
            continue
        
        # Extract key components
        split = key.split(":")
        chrom = split[0]
        start = int(split[-1].split("-")[0])
        
        # Slide window along the sequence
        for i in range(len(seq) - config.kmer_length + 1):
            # Extract 6-mer
            kmer = seq[i:i+config.kmer_length].upper()
            
            # Skip if invalid bases present
            if 'N' in kmer:
                continue
                
            # Match to get encoded value and strand
            index, strand = config.match_kmer(kmer)
            
            # Skip invalid k-mers
            if index >= max_codes:
                continue
                
            # Create result tuple
            result_tuple = (region_id, i + 1, strand, 0)
            
            # Add to categories
            categories[index].append(result_tuple)
        
        region_id += 1
    
    # Write the binary file
    with open(output_bin, 'wb') as bin_file:
        for index_id, category in enumerate(categories):
            # Write sentinel
            bin_file.write(struct.pack('I', config.encode_sentinel(index_id)))
            
            # Write sites
            for site in category:
                region, pos, strand, damage = site
                encoded = config.encode_site(region, pos, strand, damage)
                bin_file.write(struct.pack('I', encoded))
    
    # Initialize empty damage index
    damage_index = [[0, 0] for _ in range(max_codes)]
    
    # Write the damage index
    with open(output_index, 'w') as index_file:
        for index in damage_index:
            index_file.write(f"{index[0]}\t{index[1]}\n")
    
    print(f"Wrote {region_id} regions to {output_bin}")
    print(f"Created empty damage index at {output_index}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Encode 6-mers from genomic intervals")
    parser.add_argument("--genome", required=True, help="Path to the reference genome FASTA file")
    parser.add_argument("--intervals", required=True, help="Path to the intervals file")
    parser.add_argument("--regions", required=True, help="Path to the extracted regions FASTA file")
    parser.add_argument("--output-bin", required=True, help="Path to the output binary file")
    parser.add_argument("--output-index", required=True, help="Path to the output index file")
    
    args = parser.parse_args()
    
    # Load genome
    print("Loading genome...")
    genome = read_fasta(args.genome)
    
    # Load intervals
    print("Loading intervals...")
    intervals = hold_interval_indices(args.intervals)
    
    # Load regions
    print("Loading regions...")
    regions = read_fasta(args.regions)
    
    # Extract and encode 6-mers
    print("Extracting and encoding 6-mers...")
    extract_6mers(regions, intervals, genome, args.output_bin, args.output_index)
    
    print("Done!")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

def split_bed_by_score(input_file: Path, output_prefix: Path):
    """
    Split a BED file of TF binding sites into high and low scoring groups.
    
    Args:
        input_file: Path to input BED file
        output_prefix: Prefix for output files (will append _high.bed and _low.bed)
    """
    # Read the entire file
    df = pd.read_csv(input_file, sep='\t',
                     names=['chrom', 'start', 'end', 'strand', 'score'])
    
    # Get median score
    median_score = df['score'].median()
    
    # Split sites
    high_sites = df[df['score'] >= median_score]
    low_sites = df[df['score'] < median_score]
    
    # Create output files
    high_file = output_prefix.with_name(f"{output_prefix.stem}_high.bed")
    low_file = output_prefix.with_name(f"{output_prefix.stem}_low.bed")
    
    # Write output files
    high_sites.to_csv(high_file, sep='\t', header=False, index=False)
    low_sites.to_csv(low_file, sep='\t', header=False, index=False)

def process_tf(base_dir: str, tf_name: str):
    """
    Process a specific TF's binding sites.
    Expects base_dir/tf_name/tf_name_{plus,minus}.bed files.
    
    Args:
        base_dir: Path to directory containing TF folders
        tf_name: Name of the TF to process
    """
    base_path = Path(base_dir)
    tf_dir = base_path / tf_name
    
    if not tf_dir.exists():
        raise FileNotFoundError(f"TF directory not found: {tf_dir}")
    
    print(f"Processing {tf_name}...")
    
    # Process plus and minus files
    for strand in ['plus', 'minus']:
        bed_file = tf_dir / f"{tf_name}_{strand}.bed"
        if bed_file.exists():
            output_prefix = tf_dir / f"{tf_name}_{strand}"
            split_bed_by_score(bed_file, output_prefix)
        else:
            print(f"Warning: {bed_file} not found")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Split TF binding sites by MOODS score for a specific TF')
    parser.add_argument('base_dir', 
                       help='Base directory containing TF folders')
    parser.add_argument('tf_name',
                       help='Name of the TF to process (e.g., ETS_1)')
    
    args = parser.parse_args()
    process_tf(args.base_dir, args.tf_name)
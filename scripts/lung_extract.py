#!/usr/bin/env python3

import sys
import argparse

def parse_bed_line(line):
    """Parse a BED file line into its components."""
    fields = line.strip().split('\t')
    if len(fields) < 3:
        raise ValueError(f"Invalid BED line: {line}")
    
    chrom = fields[0]
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError:
        raise ValueError(f"Invalid coordinates in BED line: {line}")
    
    # Keep any additional fields
    additional_fields = fields[3:] if len(fields) > 3 else []
    
    return chrom, start, end, additional_fields

def shrink_interval(chrom, start, end, additional_fields):
    """
    Shrink an interval to ±75bp from its midpoint.
    Returns a tuple of (chrom, new_start, new_end, additional_fields).
    """
    # Calculate midpoint
    midpoint = (start + end) // 2
    
    # Calculate new coordinates
    new_start = midpoint - 75
    new_end = midpoint + 75
    
    # Ensure coordinates don't go below 0
    new_start = max(0, new_start)
    
    return chrom, new_start, new_end, additional_fields

def format_bed_line(chrom, start, end, additional_fields):
    """Format coordinates and fields into a BED file line."""
    base_fields = [chrom, str(start), str(end)]
    return '\t'.join(base_fields + additional_fields)

def process_bed_file(input_file):
    """Process each line of the BED file and print shrunken intervals."""
    try:
        for line in input_file:
            # Skip empty lines and comments
            if not line.strip() or line.startswith('#'):
                continue
                
            try:
                # Parse the line
                chrom, start, end, additional_fields = parse_bed_line(line)
                
                # Shrink the interval
                new_chrom, new_start, new_end, new_fields = shrink_interval(
                    chrom, start, end, additional_fields
                )
                
                # Print the new interval
                print(format_bed_line(new_chrom, new_start, new_end, new_fields))
                
            except ValueError as e:
                print(f"Warning: {e}", file=sys.stderr)
                continue
                
    except BrokenPipeError:
        # Handle broken pipe (e.g., when piping to head)
        sys.stderr.close()

def main():
    parser = argparse.ArgumentParser(
        description='Shrink BED intervals to ±75bp from their midpoint'
    )
    parser.add_argument(
        'input_file', 
        nargs='?', 
        type=argparse.FileType('r'), 
        default=sys.stdin,
        help='Input BED file (default: stdin)'
    )
    
    args = parser.parse_args()
    
    try:
        process_bed_file(args.input_file)
    finally:
        if args.input_file is not sys.stdin:
            args.input_file.close()

if __name__ == '__main__':
    main()
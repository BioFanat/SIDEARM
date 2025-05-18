#!/usr/bin/env python3
"""
Simplified test script for the configurable k-mer encoding/decoding system.
"""

import os
import struct
import tempfile
import argparse
from typing import Tuple, List, Dict, Any

from kmer_config import (
    KmerConfig, get_config_by_name, get_default_config,
    CPD_CONFIG, BPDE_CONFIG, MER6_CONFIG
)

def test_cpd_config():
    """Test the CPD configuration"""
    config = CPD_CONFIG
    print(f"\nTesting {config.name} configuration:")
    
    # Test k-mer matching
    test_kmers = [
        ('ACCT', 'Expected code based on position: 19'),
        ('TGGA', 'Expected code based on position: 48'),
        ('NNNN', 'Expected invalid code: 64')
    ]
    
    for kmer, expectation in test_kmers:
        code, strand = config.match_kmer(kmer)
        print(f"  {kmer} -> ({code}, {strand}) - {expectation}")
    
    # Test encoding/decoding
    test_sites = [
        (123, 456, '+', 7),
        (0, 0, '+', 0)
    ]
    
    for site in test_sites:
        region, pos, strand, damage = site
        encoded = config.encode_site(region, pos, strand, damage)
        decoded = config.decode_site(encoded)
        match = "MATCH" if decoded[1:] == (region, pos, strand, damage) else "MISMATCH"
        print(f"  Encoded ({region}, {pos}, {strand}, {damage}) -> {encoded}")
        print(f"  Decoded -> {decoded[1:]} - {match}")
    
    # Test sentinel
    sentinel = config.encode_sentinel(42)
    decoded = config.decode_site(sentinel)
    print(f"  Sentinel(42) -> {sentinel}")
    print(f"  Decoded -> {decoded}")

def test_bpde_config():
    """Test the BPDE configuration"""
    config = BPDE_CONFIG
    print(f"\nTesting {config.name} configuration:")
    
    # Test k-mer matching
    test_kmers = [
        ('AGT', 'Expected code based on position: 3'),
        ('ACT', 'Expected code after rev comp: 3'),
        ('ACA', 'Expected invalid code: 16')
    ]
    
    for kmer, expectation in test_kmers:
        code, strand = config.match_kmer(kmer)
        print(f"  {kmer} -> ({code}, {strand}) - {expectation}")
    
    # Test encoding/decoding
    test_sites = [
        (123, 456, '+', 7),
        (0, 0, '+', 0)
    ]
    
    for site in test_sites:
        region, pos, strand, damage = site
        encoded = config.encode_site(region, pos, strand, damage)
        decoded = config.decode_site(encoded)
        match = "MATCH" if decoded[1:5] == (region, pos, strand, damage) else "MISMATCH"
        print(f"  Encoded ({region}, {pos}, {strand}, {damage}) -> {encoded}")
        print(f"  Decoded -> {decoded[1:5]} - {match}")

def test_6mer_config():
    """Test the 6-mer configuration"""
    config = MER6_CONFIG
    print(f"\nTesting {config.name} configuration:")
    
    # Test k-mer matching
    test_kmers = [
        ('AAAAAA', 'Generic encoding'),
        ('NNNNNN', 'Expected invalid code')
    ]
    
    for kmer, expectation in test_kmers:
        code, strand = config.match_kmer(kmer)
        print(f"  {kmer} -> ({code}, {strand}) - {expectation}")
    
    # Test encoding/decoding with new format (arbitrary values)
    test_sites = [
        (123, 456, '+', 7),
        (0, 0, '+', 0)
    ]
    
    for site in test_sites:
        region, pos, strand, damage = site
        try:
            encoded = config.encode_site(region, pos, strand, damage)
            decoded = config.decode_site(encoded)
            print(f"  Encoded ({region}, {pos}, {strand}, {damage}) -> {encoded}")
            print(f"  Decoded -> {decoded[1:5]}")
        except Exception as e:
            print(f"  ERROR encoding/decoding ({region}, {pos}, {strand}, {damage}): {e}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Test the configurable k-mer system")
    parser.add_argument("--config", choices=["CPD", "BPDE", "6MER", "all"], default="all",
                        help="Configuration to test")
    
    args = parser.parse_args()
    
    # Run tests based on the configuration specified
    if args.config == "all" or args.config == "CPD":
        test_cpd_config()
    
    if args.config == "all" or args.config == "BPDE":
        test_bpde_config()
    
    if args.config == "all" or args.config == "6MER":
        test_6mer_config()
    
    print("\nTests completed!")

if __name__ == "__main__":
    main()
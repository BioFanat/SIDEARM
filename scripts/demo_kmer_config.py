#!/usr/bin/env python3
"""
Demonstration of the generalized k-mer configuration system.

This script shows how to use different k-mer configurations for
matching sequences and encoding/decoding site information.
"""

import argparse
from kmer_config import (
    KmerConfig, get_config_by_name, get_default_config,
    CPD_CONFIG, BPDE_CONFIG, MER6_CONFIG
)

def test_sequence_matching(seq: str):
    """
    Test matching a sequence with different k-mer configurations.
    
    Args:
        seq: DNA sequence to match
    """
    print(f"\nMatching sequence: {seq}")
    print("-" * 40)
    
    configs = [CPD_CONFIG, BPDE_CONFIG, MER6_CONFIG]
    
    for config in configs:
        print(f"{config.name} ({config.kmer_length}-mer):")
        
        # Check if the sequence has the right length
        if len(seq) != config.kmer_length:
            trimmed = seq[:config.kmer_length] if len(seq) > config.kmer_length else seq
            padded = trimmed.ljust(config.kmer_length, 'N')
            print(f"  WARNING: Sequence length mismatch (expected {config.kmer_length}, got {len(seq)})")
            print(f"  Using adjusted sequence: {padded}")
            test_seq = padded
        else:
            test_seq = seq
        
        # Match the sequence
        code, strand = config.match_kmer(test_seq)
        
        # Check if the result is valid
        is_valid = code < config._get_invalid_kmer_code()
        validity = "VALID" if is_valid else "INVALID"
        
        print(f"  Result: ({code}, {strand}) - {validity}")
        
        # For CPD, show decoded representation
        if config.name == 'CPD' and is_valid:
            try:
                decoded = config.decode_id(code)
                print(f"  Decoded: {decoded}")
            except:
                pass
        
        print()

def test_encoding_decoding():
    """Test encoding and decoding site information with different configurations."""
    print("\nEncoding/Decoding Test")
    print("-" * 40)
    
    configs = [CPD_CONFIG, BPDE_CONFIG, MER6_CONFIG]
    
    # Test values
    region = 123
    position = 456
    strand = '+'
    damage = 7
    
    print(f"Test values: region={region}, position={position}, strand={strand}, damage={damage}")
    print()
    
    for config in configs:
        print(f"{config.name} ({config.kmer_length}-mer):")
        
        # Encode
        try:
            encoded = config.encode_site(region, position, strand, damage)
            print(f"  Encoded: {encoded} (binary: {encoded:032b})")
            
            # Decode
            decoded = config.decode_site(encoded)
            print(f"  Decoded: {decoded}")
            
            # Check round-trip
            round_trip = (decoded[1], decoded[2], decoded[3], decoded[4])
            if round_trip == (region, position, strand, damage):
                print(f"  Round-trip: MATCH")
            else:
                print(f"  Round-trip: MISMATCH")
                print(f"    Expected: ({region}, {position}, {strand}, {damage})")
                print(f"    Got: {round_trip}")
        except Exception as e:
            print(f"  ERROR: {e}")
        
        print()

def test_sentinel():
    """Test sentinel encoding/decoding."""
    print("\nSentinel Test")
    print("-" * 40)
    
    configs = [CPD_CONFIG, BPDE_CONFIG, MER6_CONFIG]
    test_id = 42
    
    print(f"Test ID: {test_id}")
    print()
    
    for config in configs:
        print(f"{config.name} ({config.kmer_length}-mer):")
        
        # Encode sentinel
        sentinel = config.encode_sentinel(test_id)
        print(f"  Encoded sentinel: {sentinel} (binary: {sentinel:032b})")
        
        # Decode sentinel
        decoded = config.decode_site(sentinel)
        print(f"  Decoded sentinel: {decoded}")
        
        if decoded[0] == 'NEW' and decoded[1] == test_id:
            print(f"  Result: MATCH")
        else:
            print(f"  Result: MISMATCH")
            print(f"    Expected: ('NEW', {test_id}, ...)")
            print(f"    Got: {decoded}")
        
        print()

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Demonstration of the k-mer configuration system")
    parser.add_argument("--sequence", default="ACGTGA", help="DNA sequence to test")
    
    args = parser.parse_args()
    
    print("\nK-mer Configuration System Demonstration")
    print("=========================================")
    
    test_sequence_matching(args.sequence)
    test_encoding_decoding()
    test_sentinel()
    
    print("\nDemonstration completed!")

if __name__ == "__main__":
    main()
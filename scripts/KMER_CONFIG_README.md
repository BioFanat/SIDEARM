# Generalized k-mer Configuration for SIDEARM

This document describes the flexible k-mer encoding and decoding system implemented for SIDEARM.

## Overview

The new system allows for:

1. Supporting arbitrary k-mer lengths (3-mers, 4-mers, 5-mers, 6-mers, etc.)
2. Configuring how different damage types are encoded and decoded
3. Maintaining backward compatibility with existing binary formats
4. Extensibility for future damage types and encoding schemes

## Key Components

### KmerConfig Class

The `KmerConfig` class in `kmer_config.py` is the core of the new system. It defines:

- How k-mers are matched and encoded
- Bit allocations for different components (region, position, strand, damage)
- Encoding/decoding logic for both k-mers and binary site representations

### Predefined Configurations

The system comes with several predefined configurations:

- **CPD_CONFIG**: For CPD 4-mers with a central pyrimidine dimer (TT, TC, CT, CC)
- **BPDE_CONFIG**: For BPDE 3-mers with a central G/C
- **MER6_CONFIG**: Example configuration for 6-mers with a central dinucleotide

### Generalized Scripts

The following scripts have been updated to use the new configuration system:

- `interval_extraction_generalized.py`: Extracts k-mers from genomic intervals
- `redistribution_generalized.py`: Redistributes damage counts across k-mers
- `redistribution_mut_generalized.py`: Redistributes mutation counts across k-mers

## Using the System

### Command Line Flags

All generalized scripts accept a `--damage-type` flag to specify which configuration to use:

```bash
# For CPD 4-mers (default)
python redistribution_generalized.py --damage-type CPD [other args]

# For BPDE 3-mers
python redistribution_generalized.py --damage-type BPDE [other args]

# For custom 6-mers
python redistribution_generalized.py --damage-type 6MER [other args]
```

### Creating a New Configuration

To create a new configuration, add it to the `kmer_config.py` file:

```python
# New 5-mer configuration example
MER5_CONFIG = KmerConfig(
    name="5MER",
    kmer_length=5,
    motif_length=1,
    central_position=2,  # 0-indexed position of the central nucleotide
    max_region_bits=17,
    max_pos_bits=9,
    strand_bits=1,
    damage_bits=4,
    region_start_bit=1,
    version=1  # New format
)
```

Then add it to the `get_config_by_name` function:

```python
def get_config_by_name(name: str) -> KmerConfig:
    name = name.upper()
    if name == 'CPD':
        return CPD_CONFIG
    elif name == 'BPDE':
        return BPDE_CONFIG
    elif name == '6MER':
        return MER6_CONFIG
    elif name == '5MER':
        return MER5_CONFIG
    else:
        raise ValueError(f"No predefined configuration for {name}")
```

### Customizing k-mer Matching

For special k-mer matching logic, you can extend the `KmerConfig` class:

```python
class MyCustomConfig(KmerConfig):
    def _match_custom_kmer(self, seq: str) -> Tuple[int, str]:
        # Custom matching logic here
        # ...
        return encoded_value, strand
        
    def match_kmer(self, seq: str) -> Tuple[int, str]:
        # Override the main matching method
        return self._match_custom_kmer(seq)
```

## Bit Structure and Encoding

The binary encoding uses a 32-bit integer with the following structure by default:

- Bit 0: Sentinel bit (1 for 'NEW' markers, 0 for normal sites)
- Bits 1-17: Region ID (17 bits = 131,072 possible regions)
- Bits 18-26: Position within region (9 bits = 512 positions)
- Bit 27: Strand (0 for '+', 1 for '-')
- Bits 28-31: Damage value (4 bits = up to 15 damage units)

This structure can be customized in the `KmerConfig` initialization.

## Testing

A test script (`test_kmer_config.py`) is provided to validate the encoding/decoding system:

```bash
# Test all configurations
python test_kmer_config.py

# Test a specific configuration
python test_kmer_config.py --config CPD

# Skip integration tests
python test_kmer_config.py --skip-integration
```

## Backward Compatibility

The system maintains backward compatibility with existing binary formats by:

1. Using the same bit structure as the original implementation
2. Preserving the original sentinel marker approach
3. Using the same kmer matching logic for CPD and BPDE

However, newly created configurations (like 6MER) may use different bit layouts and won't be compatible with old parsers.

## Performance Considerations

- The encoding/decoding operations are highly optimized using bit manipulation
- For large datasets, consider using the numpy-based methods for bulk processing
- The system is designed to minimize memory usage while maintaining flexibility
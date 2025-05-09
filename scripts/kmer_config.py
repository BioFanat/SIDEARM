"""
Configurable k-mer encoding/decoding system for SIDEARM

This module provides a flexible framework for encoding genomic k-mers of different lengths
while maintaining compatibility with existing binary file formats.
"""

import numpy as np
from typing import Dict, Tuple, List, Union, Optional, Any

###############################################################################
#                            1) KMER CONFIGURATION                            #
###############################################################################

class KmerConfig:
    """
    Configuration for a specific damage type and its k-mer representation.
    
    This class defines how k-mers are encoded, decoded, and processed for a specific
    damage type (like CPD or BPDE).
    """
    
    def __init__(self, 
                 name: str,
                 kmer_length: int,
                 motif_length: int,
                 central_position: Union[int, Tuple[int, int]],
                 max_region_bits: int = 17,
                 max_pos_bits: int = 9,
                 strand_bits: int = 1,
                 damage_bits: int = 4,
                 region_start_bit: int = 1,
                 version: int = 0):
        """
        Initialize a k-mer configuration.
        
        Args:
            name: Name of the damage type (e.g., 'CPD', 'BPDE')
            kmer_length: Length of k-mer (e.g., 4 for CPD, 3 for BPDE)
            motif_length: Length of the central motif (e.g., 2 for CPD, 1 for BPDE)
            central_position: Position(s) of the motif within the k-mer (0-indexed)
                              Can be a single int or a tuple for multi-base motifs
            max_region_bits: Number of bits for region encoding (default: 17)
            max_pos_bits: Number of bits for position encoding (default: 9)
            strand_bits: Number of bits for strand encoding (default: 1)
            damage_bits: Number of bits for damage value encoding (default: 4)
            region_start_bit: First bit position for region encoding (default: 1)
            version: Version identifier for backwards compatibility (default: 0)
        """
        self.name = name
        self.kmer_length = kmer_length
        self.motif_length = motif_length
        
        # Ensure central_position is a tuple
        if isinstance(central_position, int):
            self.central_position = (central_position,)
        else:
            self.central_position = central_position
            
        # Validate central position is within k-mer
        for pos in self.central_position:
            if pos < 0 or pos >= kmer_length:
                raise ValueError(f"Central position {pos} outside of k-mer range (0-{kmer_length-1})")
        
        # Bit configuration
        self.max_region_bits = max_region_bits
        self.max_pos_bits = max_pos_bits
        self.strand_bits = strand_bits
        self.damage_bits = damage_bits
        self.region_start_bit = region_start_bit
        self.version = version
        
        # Calculate derived properties
        self.pos_start_bit = self.region_start_bit + self.max_region_bits
        self.strand_start_bit = self.pos_start_bit + self.max_pos_bits
        self.damage_start_bit = self.strand_start_bit + self.strand_bits
        self.total_bits = self.damage_start_bit + self.damage_bits
        
        # Validate total bits
        if self.total_bits > 32:
            raise ValueError(f"Total bits ({self.total_bits}) exceeds 32-bit limit")
            
        # Generate base mappings (for standard encoding)
        self._generate_base_mappings()
        
    def _generate_base_mappings(self):
        """Generate standard base-to-integer mappings for ACGT"""
        self.base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 15}
        self.int_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 15: 'N'}
        
        # For CPD, we need special mappings for the central dimer
        if self.name == 'CPD':
            self.dimer_to_int = {'CC': 0, 'CT': 1, 'TC': 2, 'TT': 3}
            self.int_to_dimer = {0: 'CC', 1: 'CT', 2: 'TC', 3: 'TT'}

    def match_kmer(self, seq: str) -> Tuple[int, str]:
        """
        Match a sequence to its encoded representation for this k-mer type.
        
        Args:
            seq: DNA sequence of length kmer_length
            
        Returns:
            Tuple of (encoded_value, strand)
            - encoded_value is an integer representing the k-mer
            - strand is '+' or '-' for the matched strand
        """
        seq = seq.upper()
        
        # Validate sequence length
        if len(seq) != self.kmer_length:
            return self._get_invalid_kmer_code(), "."
        
        # Special handlers for different damage types
        if self.name == 'CPD':
            return self._match_cpd_4mer(seq)
        elif self.name == 'BPDE':
            return self._match_bpde_3mer(seq)
        else:
            # Generic handler for other k-mer types
            return self._match_generic_kmer(seq)
            
    def _get_invalid_kmer_code(self) -> int:
        """Return the code for an invalid k-mer"""
        if self.name == 'CPD':
            return 64  # Traditional value for CPD
        elif self.name == 'BPDE':
            return 16  # Traditional value for BPDE
        else:
            # For generic k-mers, use a high value based on encoding capacity
            return 4 ** min(self.kmer_length, 3)  # Reasonable cap
    
    def _match_cpd_4mer(self, seq: str) -> Tuple[int, str]:
        """
        Create a numerical representation for each 4-mer that can be easily interpreted.

        Args:
            seq: 4-mer DNA sequence

        Returns:
            Tuple of (encoded_value, strand)
        """
        if 'T' in seq[1:3] or 'C' in seq[1:3]:
            strand = "+"
        else:
            strand = "-"
            seq = self.rev_complement(seq)

        # If any invalid bases, return error code
        if any(base not in 'ACGT' for base in seq):
            return self._get_invalid_kmer_code(), "."

        f1 = seq[0]
        dimer = seq[1:3]
        f2 = seq[3]

        # Create integer encoding in base 4 using original CPD encoding scheme
        # First base: 0=A, 1=C, 2=G, 3=T
        # Dimer: 0=CC, 1=CT, 2=TC, 3=TT
        # Last base: 0=A, 1=C, 2=G, 3=T
        try:
            # Original CPD encoding was for 4-mers like NCTN where the middle is pyrimidine dimer
            # ACCT = 1*16 + 0*4 + 3*1 = 16 + 0 + 3 = 19, but we represent it as "103" in base 4
            # which is 1*16 + 0*4 + 3*1 = 19 in base 10

            # For backward compatibility with intersection_counting.py
            if dimer in ["CC", "CT", "TC", "TT"]:
                dimer_idx = {"CC": 0, "CT": 1, "TC": 2, "TT": 3}[dimer]
                base1_idx = self.base_to_int[f1]
                base2_idx = self.base_to_int[f2]
                return base1_idx * 16 + dimer_idx * 4 + base2_idx, strand
            else:
                return self._get_invalid_kmer_code(), "."
        except KeyError:
            return self._get_invalid_kmer_code(), "."
    
    def _match_bpde_3mer(self, seq: str) -> Tuple[int, str]:
        """
        Create a numerical representation for BPDE 3-mers.

        Args:
            seq: 3-mer DNA sequence

        Returns:
            Tuple of (encoded_value, strand)
        """
        if len(seq) != 3 or seq[1] not in ['G', 'C']:
            return self._get_invalid_kmer_code(), "."

        if seq[1] == 'G':
            strand = "+"
        else:
            strand = "-"
            seq = self.rev_complement(seq)

        # If any invalid bases, return error code
        if any(base not in 'ACGT' for base in seq):
            return self._get_invalid_kmer_code(), "."

        try:
            # Original BPDE encoding for NXN where X is G
            # The 3-mer is encoded as 4*first + last where first/last are 0=A,1=C,2=G,3=T
            # For example, AGT would be 4*0 + 3 = 3
            base1_idx = self.base_to_int[seq[0]]
            base2_idx = self.base_to_int[seq[2]]
            return base1_idx * 4 + base2_idx, strand
        except KeyError:
            return self._get_invalid_kmer_code(), "."
    
    def _match_generic_kmer(self, seq: str) -> Tuple[int, str]:
        """
        Generic k-mer matcher for extensibility.
        
        For new damage types, specialize this method or add a new one.
        
        Args:
            seq: k-mer DNA sequence
            
        Returns:
            Tuple of (encoded_value, strand)
        """
        # Find the motif based on central position
        motif_pos = list(self.central_position)
        motif = ''.join(seq[pos] for pos in motif_pos)
        
        # Determine strand based on the motif
        if self._is_forward_motif(motif):
            strand = "+"
        else:
            strand = "-"
            seq = self.rev_complement(seq)
            # Need to recalculate motif after rev complement
            motif = ''.join(seq[pos] for pos in motif_pos)
        
        # Check for invalid motif
        if not self._is_valid_motif(motif):
            return self._get_invalid_kmer_code(), strand
        
        # Encode based on base composition
        encoding = 0
        base_power = 1
        for i, base in enumerate(seq):
            if base not in self.base_to_int:
                return self._get_invalid_kmer_code(), strand
            encoding += self.base_to_int[base] * base_power
            base_power *= 4
        
        return encoding, strand
    
    def _is_forward_motif(self, motif: str) -> bool:
        """
        Determine if a motif is in the forward orientation.
        Override in subclasses for specific damage types.
        """
        # Default implementation - can be overridden
        return motif in ['G', 'C', 'GG', 'CC', 'TT', 'GC', 'CG']
    
    def _is_valid_motif(self, motif: str) -> bool:
        """
        Check if a motif is valid for this damage type.
        Override in subclasses for specific damage types.
        """
        # Default implementation - can be overridden
        return all(base in 'ACGT' for base in motif)
    
    def rev_complement(self, seq: str) -> str:
        """Return the reverse complement of a DNA sequence"""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(complement.get(base, base) for base in reversed(seq))
    
    def decode_id(self, encoded_id: int) -> str:
        """
        Decode a 3-digit base-4 ID for CPDs (backward compatibility).
        
        Args:
            encoded_id: Integer ID to decode (0-63)
            
        Returns:
            Decoded sequence (e.g., "ACCT")
        """
        if self.name != 'CPD':
            raise ValueError(f"decode_id only implemented for CPD, not {self.name}")
            
        # Convert to base 4
        id_str = np.base_repr(encoded_id, base=4).zfill(3)
        
        ret = []
        ret.append(self.int_to_base.get(int(id_str[0]), 'N'))
        ret.append(self.int_to_dimer.get(int(id_str[1]), 'NN'))
        ret.append(self.int_to_base.get(int(id_str[2]), 'N'))
            
        return ''.join(ret)
    
    def encode_site(self, region: int, pos: int, strand: str, damage: int) -> int:
        """
        Encode site information into a single 32-bit integer.

        Args:
            region: Region ID (0-based)
            pos: Position within region (0-based)
            strand: Strand symbol ('+' or '-')
            damage: Damage count/value

        Returns:
            32-bit integer encoding with bit layout defined by this config
        """
        # Validate inputs
        if not (0 <= region < (1 << self.max_region_bits)):
            raise ValueError(f"Region {region} exceeds {self.max_region_bits} bits")
        if not (0 <= pos < (1 << self.max_pos_bits)):
            raise ValueError(f"Position {pos} exceeds {self.max_pos_bits} bits")
        if not (0 <= damage < (1 << self.damage_bits)):
            raise ValueError(f"Damage {damage} exceeds {self.damage_bits} bits")
        if strand not in ['+', '-']:
            raise ValueError(f"Strand must be '+' or '-', got {strand}")

        # Create binary representation
        strand_bit = 0 if strand == '+' else 1

        # For backward compatibility with original format
        if self.name == 'CPD':
            # Original CPD format used in redistribution_optimized.py:
            # bit[0] = 0 (for normal sites)
            # then skip bit[1]
            # region (17 bits)
            # position (9 bits)
            # strand (1 bit)
            # damage (3 bits)  [original used 4 bits but only 0-7 in practice]
            region_bits = f"{region:017b}"
            pos_bits = f"{pos:09b}"
            strand_bits = '0' if strand == '+' else '1'
            damage_bits = f"{damage:03b}"
            bit_str = '00' + region_bits + pos_bits + strand_bits + damage_bits + '0'
            return int(bit_str, 2)

        elif self.name == 'BPDE':
            # Original BPDE format used in redistribution_mut_optimized.py:
            # bit[0] = 0 (for normal sites)
            # region (16 bits)
            # position (9 bits)
            # damage (6 bits)
            # Note: No strand bit in original BPDE format
            region_bits = f"{region:016b}"
            pos_bits = f"{pos:09b}"
            damage_bits = f"{damage:06b}"
            bit_str = '0' + region_bits + pos_bits + damage_bits
            return int(bit_str, 2)

        # For new formats, use the configurable bit layout
        result = 0

        # Set version bit(s) if needed
        if self.version > 0:
            # For now, just set bit 0 to indicate we're using the new format
            result |= 1 << 0

        # Set region bits
        result |= (region & ((1 << self.max_region_bits) - 1)) << self.region_start_bit

        # Set position bits
        result |= (pos & ((1 << self.max_pos_bits) - 1)) << self.pos_start_bit

        # Set strand bit(s)
        result |= (strand_bit & ((1 << self.strand_bits) - 1)) << self.strand_start_bit

        # Set damage bits
        result |= (damage & ((1 << self.damage_bits) - 1)) << self.damage_start_bit

        return result
    
    def decode_site(self, encoded: int) -> Tuple[str, int, int, str, int]:
        """
        Decode a 32-bit encoded site.

        Args:
            encoded: 32-bit integer encoding

        Returns:
            Tuple of (type, region, pos, strand, damage)
            type is 'NEW' for sentinel values, 'SITE' for normal values
        """
        # Get binary representation
        binary = f'{encoded:032b}'

        # Check for sentinel values
        if int(binary[0]) == 1:
            # This is a 'NEW' sentinel
            id_value = int(binary[1:], 2)
            return ('NEW', id_value, 0, '+', 0)

        # For backward compatibility with original format
        if self.name == 'CPD':
            # Original CPD format used in redistribution_optimized.py:
            # bit[0] = 0 (for normal sites)
            # then skip bit[1]
            # region (17 bits, 2-18)
            # position (9 bits, 19-27)
            # strand (1 bit, 28)
            # damage (3 bits, 29-31)  [original used 4 bits but top bit was often unused]
            region = int(binary[2:19], 2)
            pos = int(binary[19:28], 2)
            strand_bit = int(binary[28])
            damage = int(binary[29:32], 2)
            strand = '+' if strand_bit == 0 else '-'
            return ('SITE', region, pos, strand, damage)

        elif self.name == 'BPDE':
            # Original BPDE format used in redistribution_mut_optimized.py:
            # bit[0] = 0 (for normal sites)
            # region (16 bits, 1-16)
            # position (9 bits, 17-25)
            # damage (6 bits, 26-31)
            # Note: No strand bit in original BPDE format
            region = int(binary[1:17], 2)
            pos = int(binary[17:26], 2)
            damage = int(binary[26:32], 2)
            return ('SITE', region, pos, '+', damage)  # Default to '+' strand

        # For new formats, use the configurable bit layout
        try:
            region = int(binary[self.region_start_bit:self.region_start_bit + self.max_region_bits], 2)
            pos = int(binary[self.pos_start_bit:self.pos_start_bit + self.max_pos_bits], 2)

            if self.strand_bits > 0:
                strand_bit = int(binary[self.strand_start_bit:self.strand_start_bit + self.strand_bits], 2)
                strand = '+' if strand_bit == 0 else '-'
            else:
                strand = '+'

            damage = int(binary[self.damage_start_bit:self.damage_start_bit + self.damage_bits], 2)

            return ('SITE', region, pos, strand, damage)
        except ValueError:
            # In case of parsing error, return a placeholder
            return ('ERROR', 0, 0, '+', 0)
        
    def encode_sentinel(self, id_value: int) -> int:
        """
        Create a sentinel value to mark the start of a new block.
        
        Args:
            id_value: ID for the sentinel (typically block ID)
            
        Returns:
            32-bit integer with highest bit set and ID in lower bits
        """
        # Ensure id_value fits in 31 bits
        if not (0 <= id_value < (1 << 31)):
            raise ValueError(f"ID {id_value} exceeds 31 bits")
            
        # Set highest bit and lower bits for ID
        return (1 << 31) | id_value


###############################################################################
#                       2) PREDEFINED CONFIGURATIONS                          #
###############################################################################

# CPD (4-mer) configuration
CPD_CONFIG = KmerConfig(
    name="CPD",
    kmer_length=4,
    motif_length=2,
    central_position=(1, 2),  # 0-indexed positions of the central dinucleotide
    max_region_bits=17,
    max_pos_bits=9,
    strand_bits=1,
    damage_bits=4,
    region_start_bit=1,  # Skip the first bit (0)
    version=0  # Legacy format
)

# BPDE (3-mer) configuration
BPDE_CONFIG = KmerConfig(
    name="BPDE",
    kmer_length=3,
    motif_length=1,
    central_position=1,  # 0-indexed position of the central G/C
    max_region_bits=16,
    max_pos_bits=9,
    strand_bits=0,  # Not used in original BPDE encoding
    damage_bits=6,
    region_start_bit=1,  # Skip the first bit (0)
    version=0  # Legacy format
)

# New 6-mer configuration example
MER6_CONFIG = KmerConfig(
    name="6MER",
    kmer_length=6,
    motif_length=2,
    central_position=(2, 3),  # 0-indexed positions of the central dinucleotide
    max_region_bits=17,
    max_pos_bits=9,
    strand_bits=1,
    damage_bits=4,
    region_start_bit=1,  # Skip the first bit (0)
    version=1  # New format
)


###############################################################################
#                           3) UTILITY FUNCTIONS                              #
###############################################################################

def get_config_by_name(name: str) -> KmerConfig:
    """
    Get a predefined configuration by name.
    
    Args:
        name: Name of the configuration (case-insensitive)
        
    Returns:
        KmerConfig object
    """
    name = name.upper()
    if name == 'CPD':
        return CPD_CONFIG
    elif name == 'BPDE':
        return BPDE_CONFIG
    elif name == '6MER':
        return MER6_CONFIG
    else:
        raise ValueError(f"No predefined configuration for {name}")

def get_default_config() -> KmerConfig:
    """Get the default configuration (CPD)"""
    return CPD_CONFIG


###############################################################################
#                           4) STANDALONE TESTING                             #
###############################################################################

if __name__ == "__main__":
    # Simple tests for verification
    
    # CPD test
    print("Testing CPD config:")
    config = get_config_by_name('CPD')
    seq = "ACCT"
    code, strand = config.match_kmer(seq)
    print(f"Sequence: {seq}")
    print(f"Encoded: {code}, Strand: {strand}")
    decoded = config.decode_id(code)
    print(f"Decoded: {decoded}")
    
    # Test encoding/decoding
    encoded = config.encode_site(region=123, pos=456, strand='+', damage=7)
    print(f"Encoded site: {encoded} (binary: {encoded:032b})")
    decoded = config.decode_site(encoded)
    print(f"Decoded site: {decoded}")
    
    # Test sentinel
    sentinel = config.encode_sentinel(id_value=42)
    print(f"Sentinel: {sentinel} (binary: {sentinel:032b})")
    decoded = config.decode_site(sentinel)
    print(f"Decoded sentinel: {decoded}")
    
    # BPDE test
    print("\nTesting BPDE config:")
    config = get_config_by_name('BPDE')
    seq = "AGT"
    code, strand = config.match_kmer(seq)
    print(f"Sequence: {seq}")
    print(f"Encoded: {code}, Strand: {strand}")
    
    # 6-mer test
    print("\nTesting 6-mer config:")
    config = get_config_by_name('6MER')
    seq = "ACGTGA"
    code, strand = config.match_kmer(seq)
    print(f"Sequence: {seq}")
    print(f"Encoded: {code}, Strand: {strand}")
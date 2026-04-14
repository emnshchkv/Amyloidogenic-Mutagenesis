"""
Amyloidogenic Mutagenesis Script

This script performs targeted mutagenesis to reduce amyloidogenic potential by:
1. Identifying enhancer amino acids (F, Y, W, V, L, I, Q, N, G) in specified regions
2. Replacing them with beta-breaker amino acids or sequences (R, P, WY, WM, KLVFF, LPFFD)
3. Generating combinatorial mutations and outputting results to FASTA format

"""

import argparse
import sys
import os
from itertools import combinations, product
from typing import List, Tuple, Dict, Set
import re

class AmyloidMutagenesis:
    """Main class for performing amyloidogenic mutagenesis"""
    
    ENHANCERS = set('FYWVLIQNG')
    
    BETA_BREAKERS = {
        'single': ['R', 'P'],
        'dipeptides': ['WY', 'WM'],
        'pentapeptides': ['KLVFF', 'LPFFD']
    }
    
    def __init__(
                 self, 
                 sequence: str, 
                 regions: List[Tuple[int, int]], 
                 output_file: str
                 ) -> None:
        """
        Initialize the mutagenesis object
        
        Args:
            sequence: Protein sequence string
            regions: List of tuples with (start, end) positions (1-indexed)
            output_file: Output FASTA file path
        """
        self.original_sequence = sequence.upper().strip()
        self.regions = regions
        self.output_file = output_file
        self.mutations = []
        
        # Validate sequence
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(aa in valid_aa for aa in self.original_sequence):
            raise ValueError("Invalid amino acid characters in sequence")
        
        # Validate regions
        for start, end in self.regions:
            if start < 1 or end > len(self.original_sequence) or start > end:
                raise ValueError(f"Invalid region {start}:{end} for sequence length {len(self.original_sequence)}")
    
    def find_enhancers_in_regions(self) -> Dict[int, str]:
        """
        Find enhancer amino acids in specified regions
        
        Returns:
            Dictionary mapping position (0-indexed) to amino acid
        """
        enhancers = {}
        
        for start, end in self.regions:
            # Convert to 0-indexed
            start_idx = start - 1
            end_idx = end - 1
            
            for i in range(start_idx, end_idx + 1):
                if self.original_sequence[i] in self.ENHANCERS:
                    enhancers[i] = self.original_sequence[i]
        
        return enhancers
    
    def generate_single_mutations(
                                  self, 
                                  enhancers: Dict[int, str]
                                  ) -> List[Tuple[str, str]]:
        """
        Generate single amino acid mutations
        
        Args:
            enhancers: Dictionary of enhancer positions and amino acids
            
        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        
        for pos, orig_aa in enhancers.items():
            for breaker in self.BETA_BREAKERS['single']:
                mutant_seq = list(self.original_sequence)
                mutant_seq[pos] = breaker
                mutant_sequence = ''.join(mutant_seq)
                
                description = f"{orig_aa}{pos+1}{breaker}"
                mutations.append((mutant_sequence, description))
        
        return mutations
    
    def generate_dipeptide_insertions(
                                      self, 
                                      enhancers: Dict[int, str]
                                      ) -> List[Tuple[str, str]]:
        """
        Generate dipeptide insertions at enhancer positions
        
        Args:
            enhancers: Dictionary of enhancer positions and amino acids
            
        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        
        for pos, orig_aa in enhancers.items():
            for dipeptide in self.BETA_BREAKERS['dipeptides']:
                mutant_seq = (self.original_sequence[:pos] + 
                            dipeptide + 
                            self.original_sequence[pos+1:])
                
                description = f"{orig_aa}{pos+1}_{dipeptide}"
                mutations.append((mutant_seq, description))
        
        return mutations
    
    def generate_pentapeptide_insertions(
                                         self, 
                                         enhancers: Dict[int, str]
                                         ) -> List[Tuple[str, str]]:
        """
        Generate pentapeptide insertions at enhancer positions
        
        Args:
            enhancers: Dictionary of enhancer positions and amino acids
            
        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        
        for pos, orig_aa in enhancers.items():
            for pentapeptide in self.BETA_BREAKERS['pentapeptides']:
                mutant_seq = (self.original_sequence[:pos] + 
                            pentapeptide + 
                            self.original_sequence[pos+1:])
                
                description = f"{orig_aa}{pos+1}_{pentapeptide}"
                mutations.append((mutant_seq, description))
        
        return mutations
    
    def generate_combinatorial_mutations(
                                         self, 
                                         enhancers: Dict[int, str], 
                                         max_combinations: int = 3
                                         ) -> List[Tuple[str, str]]:
        """
        Generate combinatorial mutations using single amino acid replacements
        
        Args:
            enhancers: Dictionary of enhancer positions and amino acids
            max_combinations: Maximum number of simultaneous mutations
            
        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        positions = list(enhancers.keys())
        
        for combo_size in range(2, min(len(positions) + 1, max_combinations + 1)):
            for pos_combo in combinations(positions, combo_size):
                breaker_combos = list(product(self.BETA_BREAKERS['single'], repeat=combo_size))
                
                for breakers in breaker_combos:
                    mutant_seq = list(self.original_sequence)
                    descriptions = []
                    
                    for i, (pos, breaker) in enumerate(zip(pos_combo, breakers)):
                        orig_aa = enhancers[pos]
                        mutant_seq[pos] = breaker
                        descriptions.append(f"{orig_aa}{pos+1}{breaker}")
                    
                    mutant_sequence = ''.join(mutant_seq)
                    description = "_".join(descriptions)
                    mutations.append((mutant_sequence, description))
        
        return mutations
    
    def generate_fixed_combinations(
                                    self, 
                                    enhancers: Dict[int, str]
                                    ) -> List[Tuple[str, str]]:
        """
        Generate mutations using fixed beta-breaker combinations
        
        Args:
            enhancers: Dictionary of enhancer positions and amino acids
            
        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        
        # Fixed combination 1: All enhancers to Proline
        if enhancers:
            mutant_seq = list(self.original_sequence)
            descriptions = []
            
            for pos, orig_aa in enhancers.items():
                mutant_seq[pos] = 'P'
                descriptions.append(f"{orig_aa}{pos+1}P")
            
            mutant_sequence = ''.join(mutant_seq)
            description = "ALL_TO_P_" + "_".join(descriptions)
            mutations.append((mutant_sequence, description))
        
        # Fixed combination 2: All enhancers to Arginine
        if enhancers:
            mutant_seq = list(self.original_sequence)
            descriptions = []
            
            for pos, orig_aa in enhancers.items():
                mutant_seq[pos] = 'R'
                descriptions.append(f"{orig_aa}{pos+1}R")
            
            mutant_sequence = ''.join(mutant_seq)
            description = "ALL_TO_R_" + "_".join(descriptions)
            mutations.append((mutant_sequence, description))
        
        # Fixed combination 3: Alternating P and R
        if len(enhancers) > 1:
            mutant_seq = list(self.original_sequence)
            descriptions = []
            positions = sorted(enhancers.keys())
            
            for i, pos in enumerate(positions):
                orig_aa = enhancers[pos]
                breaker = 'P' if i % 2 == 0 else 'R'
                mutant_seq[pos] = breaker
                descriptions.append(f"{orig_aa}{pos+1}{breaker}")
            
            mutant_sequence = ''.join(mutant_seq)
            description = "ALT_PR_" + "_".join(descriptions)
            mutations.append((mutant_sequence, description))
        
        return mutations
    
    def run_mutagenesis(
                        self, 
                        include_single: bool = True, 
                        include_dipeptides: bool = True, 
                        include_pentapeptides: bool = True, 
                        include_combinatorial: bool = True,
                        include_fixed: bool = True, 
                        max_combinations: int = 3
                        ) -> None:
        """
        Run the complete mutagenesis pipeline
        
        Args:
            include_single: Include single amino acid mutations
            include_dipeptides: Include dipeptide insertions
            include_pentapeptides: Include pentapeptide insertions
            include_combinatorial: Include combinatorial mutations
            include_fixed: Include fixed combination mutations
            max_combinations: Maximum number of simultaneous mutations for combinatorial
        """
        print(f"Original sequence: {self.original_sequence}")
        print(f"Original sequence length: {len(self.original_sequence)}")
        print(f"Analyzing regions: {self.regions}")
        
        enhancers = self.find_enhancers_in_regions()
        
        if not enhancers:
            print("No enhancer amino acids found in specified regions!")
            return
        
        print(f"Found {len(enhancers)} enhancer positions: {enhancers}")
        
        if include_single:
            single_muts = self.generate_single_mutations(enhancers)
            self.mutations.extend(single_muts)
            print(f"Generated {len(single_muts)} single mutations")
        
        if include_dipeptides:
            dipeptide_muts = self.generate_dipeptide_insertions(enhancers)
            self.mutations.extend(dipeptide_muts)
            print(f"Generated {len(dipeptide_muts)} dipeptide insertions")
        
        if include_pentapeptides:
            pentapeptide_muts = self.generate_pentapeptide_insertions(enhancers)
            self.mutations.extend(pentapeptide_muts)
            print(f"Generated {len(pentapeptide_muts)} pentapeptide insertions")
        
        if include_combinatorial and len(enhancers) > 1:
            combo_muts = self.generate_combinatorial_mutations(enhancers, max_combinations)
            self.mutations.extend(combo_muts)
            print(f"Generated {len(combo_muts)} combinatorial mutations")
        
        if include_fixed:
            fixed_muts = self.generate_fixed_combinations(enhancers)
            self.mutations.extend(fixed_muts)
            print(f"Generated {len(fixed_muts)} fixed combination mutations")
        
        print(f"Total mutations generated: {len(self.mutations)}")
    
    def write_fasta(self) -> None:
        """Write mutations to FASTA file"""
        try:
            with open(self.output_file, 'w') as f:
                f.write(">Original_Sequence\n")
                f.write(f"{self.original_sequence}\n\n")
                
                for i, (sequence, description) in enumerate(self.mutations, 1):
                    f.write(f">Mutant_{i:04d}_{description}\n")
                    f.write(f"{sequence}\n\n")
            
            print(f"Results written to {self.output_file}")
            
        except Exception as e:
            print(f"Error writing FASTA file: {e}")
            sys.exit(1)

#!/usr/bin/env python3

import pytest
from primer_check import check_single_primer


class TestSinglePrimerCheck:
    """Test cases for the check_single_primer function that tests multiple primer quality checks at once"""
    
    def test_good_primer_passes_all_checks(self):
        """Test a well-designed primer that should pass all quality checks"""
        # A primer with good GC content (~50%), low self-complementarity, 
        # no repeats, and proper GC clamp
        sequence = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"  # 33-mer, 51% GC, ends with CT
        result = check_single_primer(sequence)
        
        # Should return melting temperature (positive float) instead of -1
        assert isinstance(result, float)
        assert result > 0
        assert 50 < result < 80  # Reasonable Tm range
        
    def test_bad_gc_content_fails(self):
        """Test primer with GC content outside acceptable range"""
        # Too high GC content (>60%)
        sequence_high_gc = "GCGCGCGCGCGCGCGCGCGC"  # 100% GC
        result_high = check_single_primer(sequence_high_gc)
        assert result_high == -1
        
        # Too low GC content (<40%)
        sequence_low_gc = "ATATATATATATATATATAT"  # 0% GC
        result_low = check_single_primer(sequence_low_gc)
        assert result_low == -1
        
    def test_high_self_complementarity_fails(self):
        """Test primer with high self-complementarity"""
        # A palindromic sequence that will have high self-complementarity
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # Long palindrome
        result = check_single_primer(sequence)
        assert result == -1
        
    def test_bad_gc_clamp_fails(self):
        """Test primer without proper GC clamp"""
        # No GC in last 2 positions
        sequence_no_clamp = "ATCGATCGATCGATCGATAT"  # Ends with ATAT
        result_no_clamp = check_single_primer(sequence_no_clamp)
        assert result_no_clamp == -1
        
        # Too many GC in last 5 positions
        sequence_too_many_gc = "ATCGATCGATCGGCGCGCG"  # Last 5 are GCGCG (5 GCs)
        result_too_many = check_single_primer(sequence_too_many_gc)
        assert result_too_many == -1
        
    def test_repetitive_sequences_fail(self):
        """Test primers with repetitive sequences"""
        # Long single nucleotide repeat
        sequence_single_repeat = "ATCGAAAAAATCGATCGATCG"
        result_single = check_single_primer(sequence_single_repeat)
        assert result_single == -1
        
        # Dinucleotide repeat
        sequence_dinucleotide = "ATCGATATATATCGATCGATCG"
        result_dinucleotide = check_single_primer(sequence_dinucleotide)
        assert result_dinucleotide == -1
        
        # Multiple long repeats
        sequence_multiple = "ATCGTTTTTGCGCCCCTGATCG"
        result_multiple = check_single_primer(sequence_multiple)
        assert result_multiple == -1
        
    def test_case_insensitive_input(self):
        """Test that function handles lowercase input correctly"""
        sequence_upper = "CAAGAATCTGCCCACAAGGG"
        sequence_lower = "caagaatctgcccacaaggg"
        
        result_upper = check_single_primer(sequence_upper)
        result_lower = check_single_primer(sequence_lower)
        
        # Both should return the same result
        assert result_upper == result_lower
        assert isinstance(result_upper, float)
        assert result_upper > 0
        
    def test_custom_parameters(self):
        """Test function with custom parameter values"""
        sequence = "ATCGATCGATCGATCGATCG"
        
        # Test with stricter GC content requirements
        result_strict_gc = check_single_primer(sequence, min_gc=0.5, max_gc=0.5)
        # This should fail because our sequence has 50% GC but we require exactly 50%
        # (50% is not > 50% in the range check)
        assert result_strict_gc == -1
        
        # Test with more lenient self-complementarity threshold
        result_lenient_self_comp = check_single_primer(sequence, self_comp_max=20)
        assert isinstance(result_lenient_self_comp, float)
        assert result_lenient_self_comp > 0
        
    def test_edge_case_short_sequence(self):
        """Test with a shorter sequence"""
        sequence = "ACGCTCTTCCGATCT"  # 15-mer
        result = check_single_primer(sequence)
        
        # Should still work if it meets all criteria
        assert isinstance(result, float)
        assert result > 0
        
    def test_very_good_primer(self):
        """Test a very well-designed primer"""
        # A primer specifically designed to pass all checks:
        # - ~50% GC content
        # - Low self-complementarity
        # - No repeats
        # - Good GC clamp (ends with GC)
        sequence = "TGAACTTGGCCACAATCAGC"  # 20-mer, ends with GC
        result = check_single_primer(sequence)
        
        assert isinstance(result, float)
        assert result > 0
        assert 50 < result < 80  # Reasonable Tm range


if __name__ == "__main__":
    pytest.main([__file__])

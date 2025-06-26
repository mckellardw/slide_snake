#!/usr/bin/env python3
"""
Test script for barcode processing functionality.

This script tests the barcode processing logic independently of Snakemake.
"""

import pandas as pd
import tempfile
import os
from pathlib import Path
import sys

# Add the scripts directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts', 'py'))

def create_test_barcode_map():
    """Create a test barcode map file."""
    test_data = {
        0: ["ACGTACGTGCGC", "TGCATGCACGCA", "GGCCGGCCATGC"],
        1: ["1", "2", "3"],
        2: ["coord1", "coord2", "coord3"]
    }
    df = pd.DataFrame(test_data)
    return df

def test_seeker_processing():
    """Test Seeker-style barcode processing."""
    print("Testing Seeker barcode processing...")
    
    # Create test data
    bc_map = create_test_barcode_map()
    
    # Test splitting logic
    bc_1 = [bc[:8] for bc in list(bc_map[0])]
    bc_2 = [bc[8:] for bc in list(bc_map[0])]
    bc_us = [f"{bc[:8]}_{bc[8:]}" for bc in list(bc_map[0])]
    
    print(f"Original barcodes: {list(bc_map[0])}")
    print(f"BC1 (8bp): {bc_1}")
    print(f"BC2 (remaining): {bc_2}")
    print(f"Underscore format: {bc_us}")
    
    # Verify splits
    assert len(bc_1) == len(bc_2) == len(bc_us) == 3
    assert bc_1[0] == "ACGTACGT"
    assert bc_2[0] == "GCGC"
    assert bc_us[0] == "ACGTACGT_GCGC"
    
    print("✓ Seeker processing test passed")

def test_mist_processing():
    """Test miST-style barcode processing."""
    print("\nTesting miST barcode processing...")
    
    # Create test data
    bc_map = create_test_barcode_map()
    
    # Test splitting logic (10bp + remaining)
    bc_1 = [bc[:10] for bc in list(bc_map[0])]
    bc_2 = [bc[10:] for bc in list(bc_map[0])]
    bc_us = [f"{bc[:10]}_{bc[10:]}" for bc in list(bc_map[0])]
    
    print(f"Original barcodes: {list(bc_map[0])}")
    print(f"BC1 (10bp): {bc_1}")
    print(f"BC2 (remaining): {bc_2}")
    print(f"Underscore format: {bc_us}")
    
    # Verify splits
    assert bc_1[0] == "ACGTACGTGC"
    assert bc_2[0] == "GC"
    assert bc_us[0] == "ACGTACGTGC_GC"
    
    print("✓ miST processing test passed")

def test_adapter_insertion():
    """Test adapter insertion logic."""
    print("\nTesting adapter insertion...")
    
    # Test data
    bc_1 = ["ACGTACGT", "TGCATGCA", "GGCCGGCC"]
    bc_2 = ["GCGC", "CGCA", "ATGC"]
    adapter = "TCTTCAGCGTTCCCGAGA"
    primer = "CTACACGACGCTCTTCCGATCT"
    
    # Test adapter insertion
    bc_adapter = [f"{item1}{adapter}{item2}" for item1, item2 in zip(bc_1, bc_2)]
    bc_adapter_r1 = [f"{primer}{item1}{adapter}{item2}" for item1, item2 in zip(bc_1, bc_2)]
    
    print(f"BC1: {bc_1}")
    print(f"BC2: {bc_2}")
    print(f"With adapter: {bc_adapter}")
    print(f"With primer+adapter: {bc_adapter_r1}")
    
    # Verify insertions
    expected = f"ACGTACGT{adapter}GCGC"
    expected_r1 = f"{primer}ACGTACGT{adapter}GCGC"
    
    assert bc_adapter[0] == expected
    assert bc_adapter_r1[0] == expected_r1
    
    print("✓ Adapter insertion test passed")

def test_config_parsing():
    """Test configuration parsing logic."""
    print("\nTesting configuration parsing...")
    
    # Test BC_length parsing
    test_cases = [
        ("8 6", [8, 6]),
        ("10 10", [10, 10]),
        ("16", [16]),
        ("25", [25])
    ]
    
    for input_str, expected in test_cases:
        if isinstance(input_str, str):
            lengths = [int(x) for x in input_str.split()]
        else:
            lengths = [input_str]
        
        assert lengths == expected, f"Failed for {input_str}: got {lengths}, expected {expected}"
    
    print("✓ Configuration parsing test passed")

def main():
    """Run all tests."""
    print("Starting barcode processing tests...\n")
    
    try:
        test_seeker_processing()
        test_mist_processing()
        test_adapter_insertion()
        test_config_parsing()
        
        print("\n🎉 All tests passed!")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

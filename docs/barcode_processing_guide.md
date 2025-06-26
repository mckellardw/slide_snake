# Barcode Processing Configuration Guide

This document describes how the barcode processing system works and how to add support for new spatial RNA-seq technologies.

## Overview

The barcode processing system has been refactored to be more modular and extensible. Instead of hard-coded if/elif blocks for each recipe, it now uses configuration-driven processing based on the recipe sheet.

## Key Components

### 1. Recipe Sheet Configuration
The `resources/recipe_sheet.csv` file contains all the configuration needed for barcode processing:

- `BC_length`: Space-separated barcode lengths (e.g., "8 6" for two barcodes of 8bp and 6bp)
- `BC_concat`: Whether barcodes should be concatenated (True/False)
- `internal_adapter`: Adapter sequence for insertion between barcode segments
- `fwd_primer`: Forward primer sequence

### 2. Processing Scripts
- `scripts/py/process_barcode_whitelist.py`: Handles barcode splitting and whitelist generation
- `scripts/py/process_adapter_insertion.py`: Handles adapter insertion into barcodes

### 3. Helper Functions
New helper functions in `rules/0_utils.smk`:
- `get_barcode_config()`: Extracts barcode configuration from recipe sheet
- `get_recipe_barcode_strategy()`: Determines processing strategy for a recipe

## Adding New Recipes

### Step 1: Update Recipe Sheet
Add a new row to `resources/recipe_sheet.csv` with the following key fields:

```csv
recipe_name,BC_length,BC_concat,internal_adapter,fwd_primer,...
new_recipe,10 8,True,ACGTACGT,CTACACGACGCTCTTCCGATCT,...
```

### Step 2: Test Configuration
The system will automatically:
1. Parse the barcode lengths from `BC_length`
2. Determine if splitting is needed based on `BC_concat`
3. Use the appropriate processing strategy
4. Generate all required output files

### Step 3: Verify Processing
Check that the following files are generated correctly:
- `{SAMPLE}/bc/whitelist_1.txt` - First barcode segment
- `{SAMPLE}/bc/whitelist_2.txt` - Second barcode segment  
- `{SAMPLE}/bc/whitelist_underscore.txt` - Combined barcodes with underscore
- `{SAMPLE}/bc/map_underscore.txt` - Barcode map with underscore format

## Processing Strategies

### Single Barcode (BC_concat=False or single BC_length)
- Uses original barcode map as-is
- Creates simple whitelist file
- No barcode splitting

### Split Barcodes (BC_concat=True and multiple BC_length values)
- Splits barcodes at positions defined by BC_length
- Creates separate whitelists for each barcode segment
- Generates underscore-separated combined format for STAR

### Adapter Insertion (for recipes needing adapter sequences)
- Inserts internal_adapter between barcode segments
- Adds fwd_primer prefix for R1 processing
- Used by recipes like Seeker that need adapter-based alignment

## Example Configurations

### Seeker-type (Split barcodes with adapter)
```csv
BC_length: "8 6"
BC_concat: True  
internal_adapter: "TCTTCAGCGTTCCCGAGA"
fwd_primer: "CTACACGACGCTCTTCCGATCT"
```

### Visium-type (Single barcode)
```csv
BC_length: "16"
BC_concat: True
internal_adapter: ""
fwd_primer: "CTACACGACGCTCTTCCGATCT"
```

### STOmics-type (Single barcode)
```csv
BC_length: "25"
BC_concat: True
internal_adapter: ""
fwd_primer: "CTACACGACGCTCTTCCGATCT"
```

## Benefits of New System

1. **No Code Changes**: New recipes only require recipe sheet updates
2. **Consistent Processing**: All recipes use the same processing logic
3. **Easy Testing**: Can test different barcode configurations without code changes
4. **Maintainable**: Processing logic is centralized and documented
5. **Extensible**: Easy to add new processing strategies as needed

## Troubleshooting

### Common Issues
1. **Missing barcode files**: Check that BC_length and BC_concat are set correctly
2. **Wrong split positions**: Verify BC_length values match expected barcode structure
3. **Adapter insertion failures**: Ensure internal_adapter and fwd_primer are specified

### Debug Tips
1. Check the log files in `{SAMPLE}/bc/info.log` for processing details
2. Verify recipe sheet values with `get_recipe_info()` function
3. Test with a small barcode map first

## Future Enhancements

Potential improvements to consider:
1. Support for variable-length barcodes
2. Multiple adapter insertion points
3. Custom barcode validation rules
4. Integration with quality control metrics

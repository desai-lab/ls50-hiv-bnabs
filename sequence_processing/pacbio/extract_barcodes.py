#!/usr/bin/env python3
"""
extract_barcodes.py

Extracts barcodes from grab_barcodes.csv and prepares input for starcode.
Run this BEFORE submitting the starcode job.

Usage:
    python extract_barcodes.py

This will create 'barcodes.txt' with one barcode per line.
"""

import pandas as pd
import sys
from pathlib import Path

def main():
    # Input CSV file
    input_csv = '../grab_barcodes.csv'
    
    # Output file for starcode
    output_file = '../barcodes.txt'
    
    print(f"Reading {input_csv}...")
    
    if not Path(input_csv).exists():
        print(f"ERROR: File not found: {input_csv}")
        print("Make sure you're running this from the correct directory.")
        sys.exit(1)
    
    # Read the CSV and extract barcodes
    df = pd.read_csv(input_csv)
    df=df.dropna(subset='bc')
    df['bc']=df['bc'].astype(str)
    df['bc_len']=df['bc'].apply(lambda x: len(x))
    df=df[(df['bc_len']>=40) & (df['bc_len']<=70)]
    df.to_csv('../grab_barcodes_filtered.csv',index=False)
    # Get unique barcodes (drop NaN values and convert to string)
    barcodes = df['bc'].dropna().astype(str).unique()
    length_counts = (
    df['bc']
    .dropna()
    .astype(str)
    .str.len()
    .value_counts()
    .sort_index()
    .reset_index()
    )
    length_counts.columns = ['length', 'count']
    length_counts.to_csv('../length_counts.csv',index=False)
    print(f"Found {len(barcodes)} unique barcodes")
    
    # Write to output file (one per line)
    with open(output_file, 'w') as f:
        for bc in barcodes:
            f.write(f"{bc}\n")
    
    print(f"Wrote barcodes to {output_file}")
    print(f"File size: {Path(output_file).stat().st_size / 1024 / 1024:.2f} MB")
    print()
    print("Next steps:")
    print("  1. Edit run_starcode.sbatch to set INPUT_FILE='barcodes.txt'")
    print("  2. Submit with: sbatch run_starcode.sbatch")

if __name__ == "__main__":
    main()

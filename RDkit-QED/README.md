# Calculate QED Scores with calculateQED.py

## Purpose

Calculate a QED score and properties related to drug-likeness from a .sdf structure file.

## Installation

### Option 1: Install conda enviroment with provided .yml file

    conda create -f rdkit_enviroment.yml

### Option 2: Manually install dependancies

* rdkit
* openpyxl
* pandas

## Input file format

### Input files must be SDF files.
* Input files (.sdf) can be made in ChemDraw by saving as MDL SDfile V3000 (*.sdf).
* ChemDraw files with multiple structures can be saved to a single SDF file.
* All chemical structures in the input SDF files will be analyzed.

## Arguments

### Required (one of the following):

    -f (--file): Specify path of an input .sdf file.
    -d (--directory): Specify path of an input directory containing .sdf file(s).

### Optional:

    -o (--out): Specify an output directory (relative to cwd), default="QED_scores_out".

## Usage examples

### Activate rdkit conda enviroment before running

    activate rdkit

### Calculate QED score for all files in a directory

    python calculateQED.py --directory data

### Calculate QED score for all structures in a single file

    python calculateQED.py --file data/exmaple.sdf

## Output files

1. Excel file
    * File name
    * SMILES
    * QED properties
    * QED score

2. HTML file
    * File name
    * SMILES
    * QED properties
    * QED score
    * Chemical structure

![QED scores out](QED_scores_out/QED_scores_out.png)



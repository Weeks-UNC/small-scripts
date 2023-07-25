#!/usr/bin/env python3
# -----------------------------------------------------
# Calculate QED scores
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 0.1.0
#
# -----------------------------------------------------

import argparse
from glob import glob
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Chem import PandasTools
import pandas as pd


def get_QED(directory, file):
    QED_out = []

    if directory:
        sdf_files = directory.glob('*.sdf')
    if file:
        sdf_files = [file]

    for sdf in sdf_files:
        mols = [m for m in Chem.SDMolSupplier(str(sdf)) if m != None]
        for mol in mols:
            score = QED.default(mol)
            properties = QED.properties(mol)
            smi = Chem.MolToSmiles(mol)
            print()
            print(sdf, score)

            QED_out.append((sdf, smi, properties[0], properties[1], properties[2], properties[3], properties[4], properties[5], properties[6], properties[7], score))

    QED_df = pd.DataFrame(QED_out, columns=['File name', 'SMILES', 'MW', 'ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB', 'AROM', 'ALERTS', 'QED score'])
    PandasTools.AddMoleculeColumnToFrame(QED_df, smilesCol='SMILES')

    return QED_df 

def save_QED(QED_df, out):
    if not out.is_dir():
        out.mkdir(parents=True, exist_ok=True)

    QED_df.to_excel(out / 'QED_scores_out.xlsx',
            sheet_name='RDkit')
    QED_df.to_html(out / 'QED_scores_out.html')



def parseArgs():
    prs = argparse.ArgumentParser()
    file_input = prs.add_mutually_exclusive_group(required=True)
    file_input.add_argument('-d', '--directory', type=Path, required=False,
                     help='Specify a directory containing .sdf files in '
                     'order to calculate QED scores.')
    file_input.add_argument('-f', '--file', type=Path, required=False,
                     help='Specify a .sdf files in '
                     'order to calculate QED scores.')
    prs.add_argument('-o', '--out', type=Path, required=False, default=Path.cwd() / "QED_scores_out",
                     help='Specify a path for the output directory.')


    args = prs.parse_args()
    return args


def main(directory, file, out):

    QED_df = get_QED(directory, file)
    save_QED(QED_df, out)


if __name__ == "__main__":
    main(**vars(parseArgs()))

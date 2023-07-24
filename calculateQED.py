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

from glob import glob
from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Chem import PandasTools
import argparse
import pandas as pd


def get_QED(dir):
    sdfFiles = glob(f"{dir}/*.sdf")
    QED_out = []

    for sdf in sdfFiles:
        mols = [m for m in Chem.SDMolSupplier(sdf) if m != None]
        for mol in mols:
            # print(mol)
            # print(QED.default(mol))
            try:
                score = QED.default(mol)
            except:
                print(f'Skipped input file: {sdf}.\n')
                continue
            properties = QED.properties(mol)
            smi = Chem.MolToSmiles(mol)
            print(sdf, smi, score, properties)
            print()
            # print(rdkit.Chem.ALLChem.ComputeMolVolume(sdf))
            QED_out.append((sdf, smi, properties[0], properties[1], properties[2], properties[3], properties[4], properties[5], properties[6], properties[7], score))

    QED_df = pd.DataFrame(QED_out, columns=['File name', 'SMILES', 'MW', 'ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB', 'AROM', 'ALERTS', 'QED score'])
    PandasTools.AddMoleculeColumnToFrame(QED_df, smilesCol='SMILES')
    QED_df.to_excel(f"{dir}/QED_scores_out.xlsx",
             sheet_name='RDkit')
    QED_df.to_html(f"{dir}/QED_scores_out.html") 


def parseArgs():
    prs = argparse.ArgumentParser()

    prs.add_argument('-d', '--dir', type=str, required=True,
                     help='Specify a directory containing .sdf files in '
                     'order to calculate QED scores.')


    args = prs.parse_args()
    return args


def main(dir):
    get_QED(dir)


if __name__ == "__main__":
    main(**vars(parseArgs()))

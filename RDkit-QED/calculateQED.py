#!/usr/bin/env python3
# -----------------------------------------------------
# Calculate QED scores
# Seth D. Veenbaas
# Weeks Lab, UNC-CH
# 2025
#
# Version 1.2.0
#
# -----------------------------------------------------

import pandas as pd
from glob import glob
from rdkit import Chem
from rdkit.Chem import AllChem, QED, PandasTools, rdMolDescriptors, Crippen, rdDistGeom, Draw
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pickle
from PIL import Image
from PIL import PngImagePlugin

# Increase the PngImagePlugin.MAX_TEXT_CHUNK limit to 4 GB
PngImagePlugin.MAX_TEXT_CHUNK = 4 * 1024 * 1024 * 1024

### PandasTools settings
PandasTools.InstallPandasTools()
PandasTools.RenderImagesInAllDataFrames(images=True)


def make_from_smiles(
    csv: str | Path,
    out: str | Path,
    column: str,
    substructures: list | None,
    properties: bool,
    moldescriptor: bool,
    geometry: bool,
    energy_range: float,
    num_confs: int,
    no_img: bool,
):
    smiles_df = make_dataframe(csv)
    try:
        smiles_col = next(
            col for col in smiles_df.columns if column.lower() in col.lower()
        )
    except StopIteration:
        print(
            f'Error: the input file ({csv}) does not contain a column title containing "{column}"'
        )
        exit(1)
    smiles_df[smiles_col] = smiles_df[smiles_col].astype(str)

    molcol = "ROMol"
    threedcol = "3DMol"
    PandasTools.AddMoleculeColumnToFrame(
        smiles_df, smilesCol=smiles_col, molCol=molcol, includeFingerprints=True
    )

    # filters dataframe for organic molecules
    smiles_df = smiles_df.dropna(subset=[molcol])
    smiles_df = smiles_df[
        smiles_df[molcol].apply(lambda x: x.HasSubstructMatch(Chem.MolFromSmiles("C")))
    ]
    print("Generating conformers...")
    smiles_df[threedcol] = parallel_apply(smiles_df, generate_conformers_and_optimize, molcol, num_confs=num_confs)

    if substructures:
        print("Removing substructures...")
        smiles_df = remove_substructure(smiles_df, molcol, substructures=substructures)
        molcol = "Fragment"
        smiles_df["3DFragment"] = parallel_apply(smiles_df, generate_conformers_and_optimize, molcol, num_confs=num_confs)
        threedcol = "3DFragment"

    add_qed(
        smiles_df,
        molcol,
        threedcol,
        out,
        csv,
        properties,
        moldescriptor,
        geometry,
        energy_range,
        no_img,
    )


def make_from_sdf(
    sdf: str | Path,
    out: str | Path,
    substructures: list | None,
    properties: bool,
    moldescriptor: bool,
    geometry: bool,
    energy_range: float,
    num_confs: int,
    no_img: bool,
):
    name = Path(sdf).stem
    string_sdf = str(sdf)

    molcol = "ROMol"
    threedcol = "3DMol"
    sdf_df = PandasTools.LoadSDF(
        string_sdf,
        idName="ID",
        molColName=molcol,
        includeFingerprints=False,
        isomericSmiles=True,
        smilesName="SMILES",
        embedProps=True,
        removeHs=False,
        strictParsing=True,
    )

    # filters dataframe for organic molecules
    sdf_df = sdf_df.dropna(subset=[molcol])
    sdf_df = sdf_df[
        sdf_df[molcol].apply(lambda x: x.HasSubstructMatch(Chem.MolFromSmiles("C")))
    ]
    print("Generating conformers...")
    sdf_df[threedcol] = parallel_apply(sdf_df, generate_conformers_and_optimize, molcol, num_confs=num_confs)

    if substructures:
        print("Removing substructures...")
        sdf_df = remove_substructure(sdf_df, molcol, substructures=substructures)
        molcol = "Fragment"
        sdf_df["3DFragment"] = parallel_apply(sdf_df, generate_conformers_and_optimize, molcol, num_confs=num_confs)
        threedcol = "3DFragment"

    add_qed(
        sdf_df,
        molcol,
        threedcol,
        out,
        sdf,
        properties,
        moldescriptor,
        geometry,
        energy_range,
        no_img,
    )


def make_dataframe(file: Path):
    """Reads a file (CSV or XLSX) into a Pandas DataFrame.

    Args:
        file: The path to the file (.csv or .xlsx) to read.

    Returns:
        A Pandas DataFrame containing the data from the file.
    """
    extension = file.suffix.lower()

    if extension == ".csv":
        df = pd.read_csv(file)
    elif extension == ".xlsx":
        df = pd.read_excel(file)
    else:
        raise ValueError(f"Unsupported file format: {file}")

    return df


def generate_conformers_and_optimize(mol: Chem.rdchem.Mol, num_confs: int):
    try:
        mol = Chem.AddHs(mol)  # Add hydrogens
        etkdg = AllChem.ETKDGv3()  # Use ETKDG parameters
        etkdg.optimizerForceTol = float(0.0135)
        etkdg.randomSeed = 0xa700f
        etkdg.verbose = False
        etkdg.numThreads = 8
        ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=etkdg)  # Generate multiple conformers
        if not ids:  # Embedding failed
            print(f"Failed to embed molecule: {Chem.MolToSmiles(mol)}")
            return None
        for conf_id in ids:
            AllChem.UFFOptimizeMolecule(mol, confId=conf_id)  # Optimize each conformer geometry
        if mol.GetNumConformers() == 0:
            print(f"No conformers generated for molecule: {Chem.MolToSmiles(mol)}")
            return None
        return mol
    except Exception as e:
        print(f"Error generating conformers and optimizing: {e}")
        return None


def parallel_apply(df, func, column, **kwargs):
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(apply_func, df[column], [func] * len(df), [kwargs] * len(df)))
    return results


def apply_func(x, func, kwargs):
    return func(x, **kwargs)


def remove_substructure(df: pd.DataFrame, molcol: str, substructures: list | None):
    if substructures is None:
        return df
    df["Fragment"] = df[molcol]
    for substructure in substructures:
        if substructure:
            substructure = Chem.MolFromSmarts(substructure)
            df["Fragment"] = df["Fragment"].apply(lambda x: AllChem.DeleteSubstructs(x, substructure))
    df["Fragment"].apply(lambda x: AllChem.EmbedMolecule(Chem.AddHs(x)))
    df["Fragment"].apply(lambda x: Chem.SanitizeMol(x))

    return df


def add_qed(
    df: pd.DataFrame,
    molcol: str,
    threedcol: str,
    out: str | Path,
    file: str | Path,
    properties: bool,
    moldescriptor: bool,
    geometry: bool,
    energy_range: float,
    no_img: bool,
):
    df["QED"] = parallel_apply(df, QED.default, molcol)
    if properties:
        print("Calculating QED properties...")
        qed_props = parallel_apply(df, get_qed_properties, molcol)
        for prop in ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]:
            df[prop] = [props[prop] if props else None for props in qed_props]
    if moldescriptor:
        print("Calculating molecular descriptors...")
        mol_descs = parallel_apply(df, get_mol_descriptors, molcol)
        for desc in [
            "Num Ring",
            "Num Ar Ring",
            "Num ArHetcy",
            "Num Hetcy",
            "Num Hetatm",
            "Num Spiro",
            "Frac Sp3",
            "MR",
        ]:
            df[desc] = [desc_vals[desc] if desc_vals else None for desc_vals in mol_descs]
    if geometry:
        print("Calculating molecular geometries...")
        geom_props = parallel_apply(df, get_geometry_with_range, threedcol, energy_range=energy_range)
        df["NPR1"], df["NPR2"], df["Geometry"] = zip(*geom_props)

    molcol_list = ["ROMol", "3DMol"]
    if molcol not in molcol_list:
        molcol_list.append(molcol)
    if threedcol not in molcol_list:
        molcol_list.append(threedcol)

    save_df(df=df, out=out, file=file, molcol=molcol_list, no_img=no_img)


def save_df(
    df: pd.DataFrame,
    out: str | Path,
    file: str | Path,
    molcol: str | list,
    no_img: bool,
):
    if out:
        dir = Path(out)
    else:
        dir = Path(file).parent
    if not dir.is_dir():
        dir.mkdir(parents=True, exist_ok=True)
    file_name = Path(file).stem
    sdf_out = str(Path(dir, file_name + "_qed.sdf"))

    # Ensure the molecule column contains valid RDKit molecule objects
    valid_molcol = molcol[-1]
    df = df[df[valid_molcol].apply(is_valid_molecule)]

    # Write the SDF file
    PandasTools.WriteSDF(df, sdf_out, molColName=valid_molcol, properties=list(df.columns))

    print(f"Saving {Path(dir, file_name + '_qed.html')} ....")
    try:
        df.to_html(Path(dir, file_name + "_qed.html"))
    except Exception as e:
        print(f"Unable to save .HTML file: {e}")

    # Save conformations to a pickle file
    pickle_out = str(Path(dir, file_name + "_conformers.pkl"))
    with open(pickle_out, 'wb') as f:
        pickle.dump(df, f)
    print(f"Saving conformations to {pickle_out} ....")

    print(f"Saving {Path(dir, file_name + '_qed.xlsx')} ....")
    try:
        if no_img:
            df.to_excel(Path(dir, file_name + "_qed.xlsx"), index=False)
        else:
            PandasTools.SaveXlsxFromFrame(df, Path(dir, file_name + "_qed.xlsx"), molCol=molcol, size=(150, 150))
    except Exception as e:
        print(f"Unable to save .XLSX file: {e}")

    print("Done!")


def is_valid_molecule(mol):
    return mol is not None and isinstance(mol, Chem.Mol)


def get_qed_properties(mol: Chem.rdchem.Mol):
    try:
        properties = QED.properties(mol)
        return {
            "MW": properties[0],
            "ALOGP": properties[1],
            "HBA": properties[2],
            "HBD": properties[3],
            "PSA": properties[4],
            "ROTB": properties[5],
            "AROM": properties[6],
            "ALERTS": properties[7],
        }
    except Exception as e:
        print(f"Error calculating QED properties: {e}")
        return None


def get_mol_descriptors(mol: Chem.rdchem.Mol):
    try:
        return {
            "Num Ring": rdMolDescriptors.CalcNumRings(mol),
            "Num Ar Ring": rdMolDescriptors.CalcNumAromaticRings(mol),
            "Num ArHetcy": rdMolDescriptors.CalcNumAromaticHeterocycles(mol),
            "Num Hetcy": rdMolDescriptors.CalcNumHeterocycles(mol),
            "Num Hetatm": rdMolDescriptors.CalcNumHeteroatoms(mol),
            "Num Spiro": rdMolDescriptors.CalcNumSpiroAtoms(mol),
            "Frac Sp3": rdMolDescriptors.CalcFractionCSP3(mol),
            "MR": Crippen.MolMR(mol),
        }
    except Exception as e:
        print(f"Error calculating molecular descriptors: {e}")
        return None


def get_geometry(mol: Chem.rdchem.Mol, energy_range: float = 3.0):
    try:
        if mol.GetNumConformers() > 0:
            energies = []
            for conf in mol.GetConformers():
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf.GetId())
                energies.append(ff.CalcEnergy())
            min_energy = min(energies)
            boltzmann_factors = {conf.GetId(): np.exp(-(energy - min_energy) / (0.001987 * 298.15)) for conf, energy in zip(mol.GetConformers(), energies) if energy - min_energy <= energy_range}
            total_boltzmann_factor = sum(boltzmann_factors.values())
            boltzmann_weights = {conf_id: bf / total_boltzmann_factor for conf_id, bf in boltzmann_factors.items()}

            npr1_avg = sum(weight * rdMolDescriptors.CalcNPR1(mol, confId=conf_id) for conf_id, weight in boltzmann_weights.items())
            npr2_avg = sum(weight * rdMolDescriptors.CalcNPR2(mol, confId=conf_id) for conf_id, weight in boltzmann_weights.items())

            geometry = "Balanced"
            if npr1_avg - npr2_avg + 0.5 < 0:
                geometry = "Rod-like"
            elif -npr1_avg - npr2_avg + 1.5 < 0:
                geometry = "Sphere-like"
            elif npr2_avg - 0.75 < 0:
                geometry = "Disc-like"
            return npr1_avg, npr2_avg, geometry
        else:
            return None, None, None
    except Exception as e:
        print(f"Error calculating geometry: {e}")
        return None, None, None


def get_geometry_with_range(mol, energy_range):
    return get_geometry(mol, energy_range)


def parseArgs():
    prs = argparse.ArgumentParser(description="Calculate QED scores from SMILES or SDF.")
    ex_group = prs.add_mutually_exclusive_group()

    ex_group.add_argument(
        "-d", "--dir", type=Path, help="Path to a directory containing input .sdf files."
    )

    ex_group.add_argument(
        "-csv", "--csv", type=Path, help="Path to .csv or .xlsx file with a SMILES string column."
    )

    ex_group.add_argument("-sdf", "--sdf", type=Path, help="Path to an input .sdf file.")

    ex_group.add_argument("-pkl", "--pickle", type=Path, help="Path to an input .pkl file.")

    prs.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
        help="Path to output directory. (Default=directory of input file)",
    )

    prs.add_argument(
        "-c",
        "--column",
        type=str,
        default="SMILES",
        help="Name/substring of column containing SMILES strings. (Default=SMILES)",
    )

    prs.add_argument(
        "-s",
        "--substructures",
        type=str,
        nargs="*",
        default=None,
        help="Removes substructure from property calculation. \
             Use 'diazirine_handle' to remove diazarne FFF handles. \
             (Default=None)",
    )

    prs.add_argument(
        "-p",
        "--properties",
        action="store_true",
        default=False,
        help="Adds QED properties to outputs. (Default=False)",
    )

    prs.add_argument(
        "-md",
        "--moldescriptors",
        action="store_true",
        default=False,
        help="Adds QED properties to outputs. (Default=False)",
    )

    prs.add_argument(
        "-g",
        "--geometry",
        action="store_true",
        default=False,
        help="Adds NPR1, NPR2, and geometry descriptor to outputs. (Default=False)",
    )

    prs.add_argument(
        "-er",
        "--energy_range",
        type=float,
        default=3.0,
        help="Energy range in kcal/mol for Boltzmann averaging. (Default=3.0)",
    )

    prs.add_argument(
        "-cf",
        "--conformers",
        type=int,
        default=100,
        help="Number of conformers to generate. (Default=100)",
    )

    prs.add_argument(
        "-ni",
        "--no_img",
        action="store_true",
        default=False,
        help="Include 3D molecule images from the output XLSX. (Default=False)",
    )

    args = prs.parse_args()
    if args.substructures == ["diazirine_handle"]:
        args.substructures = [
            "[#8]=[#6]-[#6]-[#6]-[#6]1(-[#6]-[#6]-[#6]#[#6])-[#7]=[#7]-1",
            "[#7]-[#6]-[#6]-[#6]1(-[#6]-[#6]-[#6]#[#6])-[#7]=[#7]-1",
            "[#6]-[#6]-[#6]1(-[#6]-[#6]-[#6]#[#6])-[#7]=[#7]-1",
        ]
    return args


def main(
    dir: str | Path,
    csv: str | Path,
    sdf: str | Path,
    pickle: str | Path,
    out: str | Path,
    column: str,
    substructures: list | None,
    properties: bool,
    moldescriptors: bool,
    geometry: bool,
    energy_range: float,
    conformers: int,
    no_img: bool,
):
    if csv:
        make_from_smiles(
            csv, 
            out,
            column, 
            substructures,
            properties,
            moldescriptors,
            geometry,
            energy_range,
            conformers,
            no_img,
        )
    elif sdf:
        make_from_sdf(
            sdf, 
            out,
            substructures,
            properties,
            moldescriptors,
            geometry,
            energy_range,
            conformers,
            no_img,
        )
    elif dir:
        sdf_files = glob(f"{dir}/*.sdf")
        for sdf in sdf_files:
            make_from_sdf(
                sdf, 
                out,
                substructures,
                properties,
                moldescriptors,
                geometry,
                energy_range,
                conformers,
                no_img,
            )
    elif pickle:
        df = pd.read_pickle(pickle)
        if out:
            dir = Path(out)
        else:
            dir = Path(pickle).parent
        if not dir.is_dir():
            dir.mkdir(parents=True, exist_ok=True)
        try:
            if no_img:
                df.to_excel(Path(dir) / (pickle.stem + "_qed.xlsx"), index=False)
            else:
                PandasTools.SaveXlsxFromFrame(df, Path(dir) / (pickle.stem + "_qed.xlsx"), molCol=["ROMol", "3DMol"], size=(150, 150))
        except Exception as e:
            print(f"Unable to save .XLSX file: {e}")


if __name__ == "__main__":
    main(**vars(parseArgs()))

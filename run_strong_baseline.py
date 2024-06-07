# Runs the full strong baseline, including smina/vina docking,
# gnina rescoring, and an input conformational ensemble.
import argparse
import os
import shutil
import subprocess

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools, rdMolTransforms


def protonate_receptor_and_ligand(folder, complex_id):
    # This function is the same as in run_single_conformer_baseline.py
    if complex_id[:4] == "8CSD":
        # 8CSD is already protonated and crashes reduce, so just copy it over
        shutil.copy(
            os.path.join(folder, f"{complex_id}_protein.pdb"),
            os.path.join(folder, f"{complex_id}_H.pdb"),
        )
    else:
        protein_in = os.path.join(folder, f"{complex_id}_protein.pdb")
        protein_out = os.path.join(folder, f"{complex_id}_H.pdb")
        with open(protein_out, "w") as f:
            subprocess.run(
                ["reduce", "-BUILD", protein_in],
                stdout=f,
                stderr=subprocess.DEVNULL,
            )
    ligand_in = os.path.join(folder, f"{complex_id}_ligand_start_conf.sdf")
    ligand_out = os.path.join(folder, f"{complex_id}_ligand_start_conf_H.sdf")
    subprocess.run(["obabel", ligand_in, "-O", ligand_out, "-p", "7.4"])


def generate_conformers(folder, complex_id, num_confs=8):
    mol = Chem.MolFromMolFile(
        os.path.join(folder, f"{complex_id}_ligand_start_conf_H.sdf")
    )
    mol.RemoveAllConformers()
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=1)
    AllChem.UFFOptimizeMoleculeConfs(mol)
    with Chem.SDWriter(
        os.path.join(folder, f"{complex_id}_ligand_multiple_confs.sdf")
    ) as writer:
        for cid in range(mol.GetNumConformers()):
            writer.write(mol, confId=cid)


def run_docking(folder, complex_id):
    crystal_ligand = Chem.MolFromMolFile(
        os.path.join(folder, f"{complex_id}_ligand.sdf")
    )
    center = rdMolTransforms.ComputeCentroid(crystal_ligand.GetConformer())
    subprocess.run(
        [
            "gnina",
            "-r",
            os.path.join(folder, f"{complex_id}_H.pdb"),
            "-l",
            os.path.join(folder, f"{complex_id}_ligand_multiple_confs.sdf"),
            "-o",
            os.path.join(folder, f"{complex_id}_multiple_confs_poses.sdf"),
            "--center_x",  # bounding box matching PoseBusters methodology
            str(center.x),
            "--center_y",
            str(center.y),
            "--center_z",
            str(center.z),
            "--size_x",
            "25",
            "--size_y",
            "25",
            "--size_z",
            "25",
            "--scoring",
            "vina",
            "--exhaustiveness",
            "4",
            "--num_modes",
            "10",
            "--seed",
            "1",
        ]
    )
    # sort the poses from the multiple conformation runs, so overall best is first
    poses = PandasTools.LoadSDF(
        os.path.join(folder, f"{complex_id}_multiple_confs_poses.sdf")
    )
    poses["CNNscore"] = poses["CNNscore"].astype(float)
    gnina_order = poses.sort_values("CNNscore", ascending=False).reset_index(drop=True)
    PandasTools.WriteSDF(
        gnina_order,
        os.path.join(folder, f"{complex_id}_multiple_confs_poses.sdf"),
        properties=list(poses.columns),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base_path",
        type=str,
        help="dir containing posebusters_benchmark_set and posebusters_308_ids.csv",
    )
    args = parser.parse_args()
    base_path = args.base_path
    if base_path is None:
        base_path = os.path.dirname(os.path.abspath(__file__))
    posebusters_308 = pd.read_csv(os.path.join(base_path, "posebusters_308_ids.csv"))
    posebusters_308["complex_id"] = (
        posebusters_308["pdb_id"] + "_" + posebusters_308["ccd_id"]
    )
    for complex_id in posebusters_308["complex_id"]:
        folder = os.path.join(
            base_path, "posebusters_paper_data/posebusters_benchmark_set", complex_id
        )
        protonate_receptor_and_ligand(folder, complex_id)
        generate_conformers(folder, complex_id)
        run_docking(folder, complex_id)

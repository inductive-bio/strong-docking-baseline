# Uses the PoseBusters library to evaluate the results of the baselines.
# It saves the results in the `results` directory.
import argparse
import multiprocessing
import os
from pathlib import Path

import pandas as pd
from posebusters import PoseBusters
from rdkit import RDLogger
from rdkit.Chem import PandasTools


def bust_complex(base_path, complex_id, suffix, use_vina_order=False):
    """
    Bust predictions for the top pose of a complex.

    The use_vina_order argument can be used to re-order the poses by the
    minimizedAffinity property before busting the top pose. This is useful since
    we've saved the poses in the gnina-ranked order.
    """
    RDLogger.DisableLog("rdApp.*")
    if use_vina_order:
        # save a copy of the poses with the vina order
        poses = PandasTools.LoadSDF(
            Path(base_path, complex_id, f"{complex_id}{suffix}.sdf")
        )
        poses["minimizedAffinity"] = poses["minimizedAffinity"].astype(float)
        vina_order = poses.sort_values("minimizedAffinity").reset_index(drop=True)
        pred_file = Path(base_path, complex_id, f"{complex_id}{suffix}_vina_order.sdf")
        PandasTools.WriteSDF(
            vina_order,
            pred_file,
            properties=list(poses.columns),
        )
    else:
        pred_file = Path(base_path, complex_id, f"{complex_id}{suffix}.sdf")
    true_file = Path(base_path, complex_id, f"{complex_id}_ligand.sdf")
    receptor_file = Path(base_path, complex_id, f"{complex_id}_protein.pdb")
    buster = PoseBusters(top_n=1)
    df = buster.bust([pred_file], true_file, receptor_file)
    df = df.reset_index(drop=True).assign(complex_id=complex_id)
    return df


if __name__ == "__main__":
    RDLogger.DisableLog("rdApp.*")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base_path",
        type=str,
        help="dir containing posebusters_paper_data and posebusters_308_ids.csv",
    )
    args = parser.parse_args()
    base_path = args.base_path
    if base_path is None:
        base_path = os.path.dirname(os.path.abspath(__file__))
    posebusters_308 = pd.read_csv(os.path.join(base_path, "posebusters_308_ids.csv"))
    posebusters_308["complex_id"] = (
        posebusters_308["pdb_id"] + "_" + posebusters_308["ccd_id"]
    )
    num_processes = multiprocessing.cpu_count()
    os.makedirs(os.path.join(base_path, "results"), exist_ok=True)
    for name, args in {
        "vina_single_conformer": ("_single_conf_poses", True),
        "gnina_single_conformer": ("_single_conf_poses", False),
        "vina_multiple_conformers": ("_multiple_confs_poses", True),
        "gnina_multiple_conformers": ("_multiple_confs_poses", False),
    }.items():
        # PoseBusters is a bit slow, so parallelize
        with multiprocessing.Pool(num_processes) as p:
            results = p.starmap(
                bust_complex,
                [
                    (
                        Path(
                            base_path,
                            "posebusters_paper_data/posebusters_benchmark_set",
                        ),
                        complex_id,
                        args[0],
                        args[1],
                    )
                    for complex_id in posebusters_308["complex_id"]
                ],
            )
        full_df = pd.concat(results, ignore_index=True)
        full_df.to_csv(
            os.path.join(base_path, f"results/posebusters_results_{name}.csv"),
            index=False,
        )
        full_df = full_df.set_index("complex_id")
        print("Method:", name)
        print(
            "PB-valid:",
            full_df.all(axis=1).mean(),
            "RMSD<2:",
            full_df["rmsd_≤_2å"].mean(),
        )

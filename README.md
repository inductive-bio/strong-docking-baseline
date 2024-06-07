# Strong Docking Baseline

This repository contains the code for the baselines described in the blog post [Approaching AlphaFold 3 docking accuracy in 100 lines of code](https://www.inductive.bio/blog/strong-baseline-for-alphafold-3-docking).

The key files are:
* `run_strong_baseline.py`: This runs the full strong baseline, including smina/vina docking, gnina rescoring, and an input conformational ensemble.
* `run_single_conformer_baseline.py`: This runs smina/vina docking and gnina rescoring, using only a single starting conformation.
* `evaluate_results.py`: This script uses the PoseBusters library to evaluate the results of the baselines. It saves the results in the `results` directory.
* `posebusters_308_ids.csv`: This file includes the 308 PDB and CCD IDs used in the final [PoseBusters paper](https://arxiv.org/pdf/2308.05777) and defined in supplementary section S5. In the AlphaFold 3 paper, this is referred to as "PoseBusters v2". "PoseBusters v1" includes an additional 120 complexes that were excluded from the final benchmark due to having crystal contacts (see [note at top of the PoseBusters Zenodo record](https://zenodo.org/records/8278563)).

## Instructions to run

### Data

Retrieve the posebusters benchmark set from [here](https://zenodo.org/records/8278563), e.g. via
```
# navigate to the base directory of this repository before running these commands
wget https://zenodo.org/records/8278563/files/posebusters_paper_data.zip
unzip posebusters_paper_data.zip -d posebusters_paper_data
```

### Docker setup

To ensure reproducibility, we've provided a docker image with all needed dependencies.

You can retrieve the docker image from Docker Hub with the following command:

```
docker pull inductivebio/docking-baseline:1.0.0
```

If you'll be running on a GPU-enabled machine (recommended), ensure you've installed the [nvidia container toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) so that the docker container can access your GPU. If you encounter issues with the GPU disconnecting, try the solution [here](https://github.com/NVIDIA/nvidia-container-toolkit/issues/48).

### Non-docker setup

If you prefer not to use Docker, you can run the code directly on your machine.
Ensure you have [gnina v1.1](https://github.com/gnina/gnina/releases/tag/v1.1),
[OpenBabel v3.1.1](https://openbabel.org/docs/Installation/install.html) and [reduce v4.14](https://github.com/rlabduke/reduce)
available in your PATH.

Install dependencies from requirements.txt into a virtual environment:

```
pip install --upgrade pip && pip install -r requirements.txt
```

The code was tested with Python 3.10 on Ubuntu 22.04.

### Run the code

The code can be run from the base directory of this repo with the following commands. Note that it will take several hours to execute (and longer without a gpu). Pose output files will be added into complex directories within `posebusters_paper_data`, and final results will be saved in the `results` directory.

With docker:
```
# generate results for the strong baseline (gnina + conformational ensemble)
# as well as vina + conformational ensemble
docker run --gpus all -v $(pwd):/src inductivebio/docking-baseline:1.0.0 python3 /src/run_strong_baseline.py


# generate results for gnina baseline (no conformational ensemble)
# as well as vina baseline
docker run --gpus all -v $(pwd):/src inductivebio/docking-baseline:1.0.0 python3 /src/run_single_conformer_baseline.py

# use posebusters to evaluate the results, saving outputs in ./results
docker run --gpus all -v $(pwd):/src inductivebio/docking-baseline:1.0.0 python3 /src/evaluate_results.py
```

Without docker:

```
# generate results for the strong baseline (gnina + conformational ensemble)
# as well as vina + conformational ensemble
python run_strong_baseline.py

# generate results for gnina baseline (no conformational ensemble)
# as well as vina baseline
python run_single_conformer_baseline.py

# use posebusters to evaluate the results, saving outputs in ./results
python evaluate_results.py
```

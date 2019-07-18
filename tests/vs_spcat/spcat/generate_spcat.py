
import os
import numpy as np
import tempfile
import shutil
import json
from pathlib import Path
from subprocess import run, PIPE


def run_spcat(name, temp_dict, constants_dict):
    """
    Function to run SPCAT from within Python. The function does so with
    temporary folders within a context manager, such that there are no
    scratch files left behind when the process is finished.
    
    Parameters
    ----------
    name : str
        Name of the molecule to use for the LineList creation
    temp_dict : dict
        Dictionary containing 
    """
    curdir = os.getcwd()
    with tempfile.TemporaryDirectory() as temp_dir:
        os.chdir(temp_dir)
        try:
            with open("molecule.var", "w+") as write_file:
                write_file.write(
                    temp_dict[".var"].format_map(constants_dict)
                )
            with open("molecule.int", "w+") as write_file:
                write_file.write(
                    temp_dict[".int"].format_map(constants_dict)
                )
        except KeyError:
            raise Exception("temp_dict is missing entries!")
        # Call SPCAT
        _ = run(["spcat", "molecule"], stdout=PIPE)
        shutil.copy("molecule.cat", Path(curdir).joinpath(f"{name}.cat"))
    os.chdir(curdir)

def generate_constants():
    C, B, A = np.sort(np.random.uniform(1000., 10000., 3))
    # Generate dipole moment as boolean values
    dipoles = np.ones(3)
    u_a, u_b, u_c = dipoles
    constants_dict = {
        "A": A,
        "B": B,
        "C": C,
        "u_a": u_a,
        "u_b": u_b,
        "u_c": u_c
    }
    return constants_dict


def read_templates(temp_path):
    temp_path = Path(temp_path)
    temp_dict = dict()
    for path in temp_path.rglob("template*"):
        temp_dict[path.suffix] = path.read_text()
    return temp_dict


def iteration(index, temp_dict):
    constants_dict = generate_constants()
    run_spcat(f"S{index:d}", temp_dict, constants_dict)
    with open(f"S{index:d}.json", "w+") as write_file:
        json.dump(constants_dict, write_file)


def main(niter=5, temp_path="../data/pickett/synthetic", out_path="outputs/dataset.pkl"):
    np.random.seed()
    temp_dict = read_templates(temp_path)
    # Set up an iterator for the batch
    iterator = range(niter)
    _ = [iteration(index, temp_dict) for index in iterator]

if __name__ == "__main__":
    results = main(200, ".")


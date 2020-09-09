
#!/bin/usr/env python

import numpy as np
import pandas as pd
import json
from pyfitter.pyfitter import AsymmetricMolecule
from pathlib import Path


def main():
    ETP =  "/Users/pcarroll/Documents/GitHub/Fitter/scripts/EigenVals_P10/J0_25_dk5.dat"
    session = AsymmetricMolecule(
            **{"et_path": ETP}
        )
    errors = list()
    counter = 0
    for path in Path().rglob("./spcat/S*.json"):
        data = json.loads(path.read_text())
        constants = [data["A"], data["B"], data["C"]]
        print(f"Simulating {constants}")
        session.simulate(constants)
        print("Done.")
        brandon_df = session.format_results()
        spcat_df = read_spcat(f"spcat/{path.stem}.cat")
        tests, error, worst_trans = compare_datasets(
            brandon_df, 
            spcat_df,
            thres=1e-4,
            ntests=50
            )
        errors.append(error)
        counter += 1
        print(f"Largest error for set: {error:.4f}")
        print(f"{worst_trans}")
    avg_error = np.average(errors)
    std_error = np.std(errors)
    print ("================")
    print ("")
    print (ETP)
    print(f"Aggregate error for {counter} tests")
    print(f"{avg_error:.6e} +/- {std_error:.6e} MHz")


def read_spcat(filepath):
    df = pd.read_fwf(
        filepath, 
        widths=[13,8,8,2,10,3,7,4,2,2,2,8,2,2], 
        header=None
        )
    columns = [
        "Frequency",
        "Uncertainty",
        "Intensity",
        "DoF",
        "Lower state energy",
        "Degeneracy",
        "ID",
        "Coding",
        "J'",
        "Ka'",
        "Kc'",
        "J''",
        "Ka''",
        "Kc''"]
    df.columns = columns
    return df


def compare_datasets(brandon_df, spcat_df, thres=1e-4, ntests=50):
    brandon_df = brandon_df.loc[brandon_df["J'"] >= 6]
    brandon_slice = brandon_df.sample(ntests, replace=False)
    tests = list()
    largest_error = 0.
    for index, row in brandon_slice.iterrows():
        spcat_slice = spcat_df.loc[
            (spcat_df["J'"] == row["J'"]) &
            (spcat_df["Ka'"] == row["Ka'"]) &
            (spcat_df["Kc'"] == row["Kc'"]) &
            (spcat_df["J''"] == row["J''"]) &
            (spcat_df["Ka''"] == row["Ka''"]) &
            (spcat_df["Kc''"] == row["Kc''"])
            ]
        if len(spcat_slice) != 0:
            error = np.abs(spcat_slice["Frequency"].values - row["frequency"])
            test_data = list(row[["J'", "Ka'", "Kc'", "J''", "Ka''", "Kc''"]])
            test_data.extend([error, error <= thres])
            if error[0] > largest_error:
                largest_error = error[0]
                worst_trans = test_data
    return tests, largest_error, worst_trans


if __name__ == "__main__":
    main()

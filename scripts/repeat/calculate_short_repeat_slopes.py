#!/usr/bin/env python3
"""
short_repeat_slope_regression.py

Estimate regression slopes between short repeat values and log10-transformed
E-values, and compute 95% confidence intervals and P-values for each species.

Input:
    short_repeat_data.xlsx

Output:
    short_repeat_slopes_CI.xlsx

Method:
    Ordinary least squares (OLS) regression:
        y ~ log10(E-value)

Author:
    Your Name
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm


# ======================
# File paths
# ======================
INPUT_FILE = "short_repeat_data.xlsx"
OUTPUT_FILE = "short_repeat_slopes_CI.xlsx"


def main():

    # ----------------------
    # Load data
    # ----------------------
    df = pd.read_excel(INPUT_FILE)

    species_col = df.columns[0]

    # ----------------------
    # Identify numeric E-value columns from header
    # ----------------------
    evalue_cols = []
    E_vals = []

    for c in df.columns[1:]:
        try:
            val = float(c)
            evalue_cols.append(c)
            E_vals.append(val)
        except:
            continue

    E_vals = np.array(E_vals)
    logE = np.log10(E_vals)

    # ----------------------
    # Regression analysis
    # ----------------------
    results = []

    for _, row in df.iterrows():

        species = row[species_col]
        y = row[evalue_cols].astype(float).values

        X = sm.add_constant(logE)
        model = sm.OLS(y, X).fit()

        slope = model.params[1]
        se = model.bse[1]
        ci_lower = slope - 1.96 * se
        ci_upper = slope + 1.96 * se
        pval = model.pvalues[1]

        results.append([
            species,
            slope,
            ci_lower,
            ci_upper,
            pval
        ])

    # ----------------------
    # Save results
    # ----------------------
    df_results = pd.DataFrame(
        results,
        columns=[
            "Species",
            "Slope",
            "CI_lower",
            "CI_upper",
            "P_value"
        ]
    )

    df_results.to_excel(OUTPUT_FILE, index=False)

    print("Analysis finished successfully.")
    print("Results saved to:", OUTPUT_FILE)
    print(df_results)


if __name__ == "__main__":
    main()

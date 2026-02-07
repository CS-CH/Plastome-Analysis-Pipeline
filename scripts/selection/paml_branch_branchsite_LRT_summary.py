#!/usr/bin/env python3
import os
import glob
import re
import csv
from scipy.stats import chi2

OUT_DIR = "results_model/"
CSV_FILE = os.path.join(OUT_DIR, "branch_branchsite_LRT.csv")
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------- 函数 -----------------------------
def extract_lnL_and_np(out_file):
    lnL, np_val = "NA", "NA"
    with open(out_file) as f:
        for line in f:
            line = line.strip()
            m_lnL = re.search(r'lnL.*:\s*(-?\d+\.\d+)', line)
            if m_lnL:
                lnL = float(m_lnL.group(1))
            m_np = re.search(r'np:\s*(\d+)', line)
            if m_np:
                np_val = int(m_np.group(1))
    return lnL, np_val

def extract_omega_M0(out_file):
    with open(out_file) as f:
        for line in f:
            if "omega (dN/dS)" in line:
                m = re.search(r"= *([\d\.]+)", line)
                if m:
                    return float(m.group(1))
    return "NA"

def extract_omega_branch(out_file):
    omega_list = []
    with open(out_file) as f:
        capture = False
        for line in f:
            if 'w (dN/dS) for branches:' in line:
                capture = True
                parts = line.split(":",1)[1].split()
                omega_list.extend([float(x) for x in parts])
            elif capture and line.strip() == "":
                break
    if len(omega_list) == 2:
        return omega_list
    return ["NA","NA"]

def extract_protein_length(out_file):
    with open(out_file) as f:
        content = f.read()
        m = re.search(r'ls\s*=\s*(\d+)', content)
        if m:
            return int(m.group(1))
    return "NA"

def calc_LRT(lnL1, lnL2, np1, np2):
    if "NA" not in [lnL1, lnL2, np1, np2]:
        df = np1 - np2
        LRT = 2 * (lnL1 - lnL2)
        p = 1 - chi2.cdf(LRT, df)
        return round(LRT,6), p
    return "NA", "NA"

# ----------------------------- 主程序 -----------------------------
genes = [d for d in glob.glob(os.path.join(OUT_DIR, "*")) if os.path.isdir(d)]

with open(CSV_FILE, "w", newline="") as f_csv:
    writer = csv.writer(f_csv)
    writer.writerow([
        "Gene",
        "lnL_M0","np_M0","omega_M0",
        "lnL_b_neut","np_b_neut","omega_b_neut",
        "lnL_b_free","np_b_free","omega_b_free_bg","omega_b_free_fg",
        "LRT_bfree_vs_M0","p_bfree_vs_M0",
        "LRT_bfree_vs_b_neut","p_bfree_vs_b_neut",
        "lnL_bsA1","np_bsA1",
        "lnL_bsA","np_bsA",
        "LRT_bsA_vs_bsA1","p_bsA_vs_bsA1",
        "Constraint_status",
        "Delta_omega_fg_bg",
        "Constraint_relaxation",
        "Branchsite_positive_selection",
        "BEB_sites",
        "BEB_sites_ratio",
        "Protein_length",
        "Abnormal_omega_flag"
    ])

    for gene_dir in genes:
        gene = os.path.basename(gene_dir)

        paths = {m: os.path.join(gene_dir, m, "out")
                 for m in ["M0","b_neut","b_free","bsA1","bsA"]}
        paths = {k:v if os.path.exists(v) else None for k,v in paths.items()}

        lnL, np_model, omega = {}, {}, {}

        for model, path in paths.items():
            if path:
                lnL[model], np_model[model] = extract_lnL_and_np(path)
                if model == "M0":
                    omega[model] = extract_omega_M0(path)
                elif model == "b_neut":
                    omega[model] = extract_omega_branch(path)[0]
                elif model == "b_free":
                    omega[model] = extract_omega_branch(path)
            else:
                lnL[model], np_model[model] = "NA", "NA"
                omega[model] = ["NA","NA"] if model=="b_free" else "NA"

        bg_omega, fg_omega = omega["b_free"]

        LRT_bfree_vs_M0, p_bfree_vs_M0 = calc_LRT(
            lnL["b_free"], lnL["M0"], np_model["b_free"], np_model["M0"])

        LRT_bfree_vs_b_neut, p_bfree_vs_b_neut = calc_LRT(
            lnL["b_free"], lnL["b_neut"], np_model["b_free"], np_model["b_neut"])

        LRT_bsA_vs_bsA1, p_bsA_vs_bsA1 = calc_LRT(
            lnL["bsA"], lnL["bsA1"], np_model["bsA"], np_model["bsA1"])

        # Constraint status
        if p_bfree_vs_M0 != "NA" and p_bfree_vs_M0 < 0.05:
            if fg_omega > bg_omega:
                Constraint_status = "Relaxed_constraint"
            elif fg_omega < bg_omega:
                Constraint_status = "Strengthened_constraint"
            else:
                Constraint_status = "No_change"
        else:
            Constraint_status = "No_change"

        # Delta omega
        Delta_omega_fg_bg = "NA"
        if Constraint_status != "No_change" and bg_omega not in ["NA",0]:
            Delta_omega_fg_bg = round((fg_omega - bg_omega)/bg_omega, 6)

        Constraint_relaxation = "Yes" if (
            p_bfree_vs_M0 != "NA" and p_bfree_vs_M0 < 0.05 and
            p_bfree_vs_b_neut != "NA" and p_bfree_vs_b_neut >= 0.05
        ) else "No"

        Branchsite_positive_selection = "Yes" if (
            p_bsA_vs_bsA1 != "NA" and p_bsA_vs_bsA1 < 0.05
        ) else "No"

        # BEB
        BEB_sites = []
        protein_length = extract_protein_length(paths["M0"]) if paths["M0"] else "NA"
        BEB_sites_ratio = ""

        if Branchsite_positive_selection == "Yes" and paths["bsA"]:
            with open(paths["bsA"]) as f:
                capture = False
                for line in f:
                    if "Bayes Empirical Bayes" in line:
                        capture = True
                    elif capture and re.match(r'^\s*\d+', line):
                        parts = line.split()
                        if parts[2].endswith("*") or float(parts[2]) >= 0.95:
                            BEB_sites.append(parts[0] + parts[1])
                    elif capture and line.strip() == "":
                        break
            if protein_length != "NA" and protein_length > 0:
                BEB_sites_ratio = round(len(BEB_sites)/protein_length, 4)

        # Abnormal omega
        Abnormal_omega_flag = "Yes" if any(
            x not in ["NA"] and (x in [0.001,999] or x > 10)
            for x in [omega["M0"], omega["b_neut"], bg_omega, fg_omega]
        ) else "No"

        writer.writerow([
            gene,
            lnL["M0"], np_model["M0"], omega["M0"],
            lnL["b_neut"], np_model["b_neut"], omega["b_neut"],
            lnL["b_free"], np_model["b_free"], bg_omega, fg_omega,
            LRT_bfree_vs_M0, p_bfree_vs_M0,
            LRT_bfree_vs_b_neut, p_bfree_vs_b_neut,
            lnL["bsA1"], np_model["bsA1"],
            lnL["bsA"], np_model["bsA"],
            LRT_bsA_vs_bsA1, p_bsA_vs_bsA1,
            Constraint_status,
            Delta_omega_fg_bg,
            Constraint_relaxation,
            Branchsite_positive_selection,
            ";".join(BEB_sites),
            BEB_sites_ratio,
            protein_length,
            Abnormal_omega_flag
        ])

        print(f"{gene}: done")

print(f"\nAll results saved to {CSV_FILE}")


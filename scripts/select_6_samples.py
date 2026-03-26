#!/usr/bin/env python3
import re
import sys
import pandas as pd

GEO_MATRIX = "meta/GSE95587_series_matrix.txt"
ENA_TSV    = "meta/PRJNA377568_ena_read_run.tsv"
OUT_TSV    = "meta/selected_6_runs.tsv"
OUT_TXT    = "meta/selected_6_runs.txt"

def die(msg: str, code: int = 1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def find_field(char_list, patterns):
    for c in char_list:
        for pat in patterns:
            if re.search(pat, c, re.I):
                return c.split(":", 1)[1].strip() if ":" in c else c.strip()
    return ""

def parse_age(x):
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else None

def norm_dx(s):
    s = str(s).lower()
    if "alzheimer" in s or s.strip() == "ad":
        return "AD"
    if "control" in s or "normal" in s or "neurologically normal" in s:
        return "Control"
    return "Other"

def main():
    # ---- 1) Parse GSM -> diagnosis, age from GEO series matrix ----
    try:
        lines = open(GEO_MATRIX, encoding="utf-8", errors="ignore").read().splitlines()
    except FileNotFoundError:
        die(f"Missing {GEO_MATRIX}. Download and gunzip the series matrix first.")

    gsm_line = next((l for l in lines if l.startswith("!Sample_geo_accession")), None)
    if gsm_line is None:
        die("Could not find !Sample_geo_accession line in GEO matrix.")
    gsms = gsm_line.split("\t")[1:]

    char_lines = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    if not char_lines:
        die("Could not find any !Sample_characteristics_ch1 lines in GEO matrix.")

    chars = [[] for _ in gsms]
    for l in char_lines:
        vals = l.split("\t")[1:]
        for i, v in enumerate(vals):
            chars[i].append(v.strip())

    rows=[]
    for gsm, clist in zip(gsms, chars):
        diagnosis_raw = find_field(clist, [r"diagnos", r"disease", r"phenotype", r"group"])
        age_raw       = find_field(clist, [r"\bage\b", r"age \(years\)", r"donor age", r"age_years"])
        rows.append((gsm, diagnosis_raw, age_raw))

    dx = pd.DataFrame(rows, columns=["GSM","diagnosis_raw","age_raw"])
    dx["age"] = dx["age_raw"].map(parse_age)
    dx["diagnosis"] = dx["diagnosis_raw"].map(norm_dx)
    dx = dx[dx["diagnosis"].isin(["AD","Control"])].dropna(subset=["age"])

    if dx.empty:
        die("No AD/Control samples with parsable age found in GEO matrix. Check matrix format/fields.")

    # ---- 2) Read ENA run table ----
    try:
        ena = pd.read_csv(ENA_TSV, sep="\t")
    except FileNotFoundError:
        die(f"Missing {ENA_TSV}. Download ENA filereport first.")

    required = {"run_accession", "fastq_ftp"}
    if not required.issubset(set(ena.columns)):
        die(f"ENA table missing required columns {required}. Found: {list(ena.columns)}")

    # columns that might contain GSM IDs
    possible_gsm_cols = [c for c in ena.columns if c.lower() in [
        "sample_alias", "secondary_sample_accession", "experiment_alias"
    ]]

    if not possible_gsm_cols:
        die(
            "ENA table doesn't include sample_alias/secondary_sample_accession/experiment_alias fields.\n"
            "Re-download ENA file with those fields included."
        )

    # build run -> (ftp, GSM) by searching GSM\d+ in those columns
    run_rows=[]
    for _, r in ena.iterrows():
        run = r.get("run_accession")
        ftp = str(r.get("fastq_ftp", ""))
        gsm_found = None
        for c in possible_gsm_cols:
            v = str(r.get(c, ""))
            m = re.search(r"(GSM\d+)", v)
            if m:
                gsm_found = m.group(1)
                break
        run_rows.append((run, ftp, gsm_found))

    runmap = pd.DataFrame(run_rows, columns=["run_accession","fastq_ftp","GSM"])
    runmap = runmap.dropna(subset=["run_accession"])
    # Keep only rows where we actually found a GSM
    runmap = runmap.dropna(subset=["GSM"])

    if runmap.empty:
        die(
            "Could not find any GSM IDs inside ENA fields (sample_alias/secondary_sample_accession/experiment_alias).\n"
            "Open meta/PRJNA377568_ena_read_run.tsv and check which field contains GSM, then request that field in ENA API."
        )

    # ---- 3) Join dx (GSM) with runmap (SRR/ENA runs) ----
    j = dx.merge(runmap, on="GSM", how="inner")
    if j.empty:
        die(
            "Join produced 0 rows (GSMs from GEO did not match GSMs parsed from ENA fields).\n"
            "Inspect GEO GSM list and ENA file; you may need a different ENA field."
        )

    # ---- 4) Select 3 AD across age range + nearest-age controls ----
    ad = j[j.diagnosis=="AD"].sort_values("age").reset_index(drop=True)
    ct = j[j.diagnosis=="Control"].sort_values("age").reset_index(drop=True)

    if len(ad) < 3 or len(ct) < 3:
        die(f"Not enough runs after join. AD runs={len(ad)}, Control runs={len(ct)}")

    idx = [0, len(ad)//2, len(ad)-1]
    ad_pick = ad.iloc[idx].copy()

    used=set()
    ctl_rows=[]
    for _, r in ad_pick.iterrows():
        ct2 = ct[~ct.run_accession.isin(used)].copy()
        ct2["agediff"] = (ct2["age"] - r["age"]).abs()
        best = ct2.sort_values(["agediff","age"]).iloc[0]
        used.add(best.run_accession)
        ctl_rows.append(best)

    ctl_pick = pd.DataFrame(ctl_rows)

    sel = pd.concat([
        ad_pick.assign(group="AD_pick"),
        ctl_pick.assign(group="Control_pick")
    ], ignore_index=True)

    sel = sel[["group","run_accession","fastq_ftp","GSM","diagnosis","age"]].sort_values(["group","age"])
    sel.to_csv(OUT_TSV, sep="\t", index=False)
    sel["run_accession"].to_csv(OUT_TXT, index=False, header=False)

    print(sel.to_string(index=False))
    print(f"\nWrote {OUT_TSV} and {OUT_TXT}")

if __name__ == "__main__":
    main()

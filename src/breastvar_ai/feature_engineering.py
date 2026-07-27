"""
Feature Engineering Module for BreastVarAI
Computes physicochemical amino acid shifts (Grantham distance, charge delta, 
polarity shift, molecular weight change) for parsed missense variants.
"""

import os
import glob
import pandas as pd

# Physicochemical Amino Acid Properties
AA_PROPERTIES = {
    "A": {"charge": 0, "polarity": "nonpolar", "weight": 89.09},
    "R": {"charge": 1, "polarity": "basic", "weight": 174.20},
    "N": {"charge": 0, "polarity": "polar", "weight": 132.12},
    "D": {"charge": -1, "polarity": "acidic", "weight": 133.10},
    "C": {"charge": 0, "polarity": "polar", "weight": 121.16},
    "E": {"charge": -1, "polarity": "acidic", "weight": 147.13},
    "Q": {"charge": 0, "polarity": "polar", "weight": 146.15},
    "G": {"charge": 0, "polarity": "nonpolar", "weight": 75.07},
    "H": {"charge": 1, "polarity": "basic", "weight": 155.16},
    "I": {"charge": 0, "polarity": "nonpolar", "weight": 131.17},
    "L": {"charge": 0, "polarity": "nonpolar", "weight": 131.17},
    "K": {"charge": 1, "polarity": "basic", "weight": 146.19},
    "M": {"charge": 0, "polarity": "nonpolar", "weight": 149.21},
    "F": {"charge": 0, "polarity": "nonpolar", "weight": 165.19},
    "P": {"charge": 0, "polarity": "nonpolar", "weight": 115.13},
    "S": {"charge": 0, "polarity": "polar", "weight": 105.09},
    "T": {"charge": 0, "polarity": "polar", "weight": 119.12},
    "W": {"charge": 0, "polarity": "nonpolar", "weight": 204.23},
    "Y": {"charge": 0, "polarity": "polar", "weight": 181.19},
    "V": {"charge": 0, "polarity": "nonpolar", "weight": 117.15}
}

# Grantham Distance Matrix Subset
GRANTHAM_MATRIX = {
    ("R", "H"): 29, ("R", "Q"): 43, ("V", "A"): 64, ("A", "E"): 107,
    ("G", "S"): 56, ("L", "P"): 98, ("R", "C"): 180, ("D", "N"): 23
}


def get_grantham_score(ref: str, alt: str) -> int:
    """Returns Grantham score for substitution, checking both direction pairs."""
    return GRANTHAM_MATRIX.get((ref, alt), GRANTHAM_MATRIX.get((alt, ref), 75))


def engineer_features(df: pd.DataFrame) -> pd.DataFrame:
    """Calculates charge delta, weight change, polarity shift, and Grantham score."""
    if df.empty:
        return df

    engineered_rows = []
    for _, row in df.iterrows():
        ref = row["ref_aa"]
        alt = row["alt_aa"]

        ref_prop = AA_PROPERTIES.get(ref, {})
        alt_prop = AA_PROPERTIES.get(alt, {})

        charge_delta = abs(alt_prop.get("charge", 0) - ref_prop.get("charge", 0))
        weight_delta = abs(alt_prop.get("weight", 0.0) - ref_prop.get("weight", 0.0))
        polarity_changed = int(ref_prop.get("polarity") != alt_prop.get("polarity"))
        grantham = get_grantham_score(ref, alt)

        target = 1 if "Pathogenic" in str(row["clinical_significance"]) else 0

        engineered_rows.append(
            {
                "clinvar_id": row["clinvar_id"],
                "gene_symbol": row["gene_symbol"],
                "protein_position": row["protein_position"],
                "ref_aa": ref,
                "alt_aa": alt,
                "charge_delta": charge_delta,
                "weight_delta": round(weight_delta, 2),
                "polarity_changed": polarity_changed,
                "grantham_score": grantham,
                "clinical_significance": row["clinical_significance"],
                "target_class": target,
            }
        )

    return pd.DataFrame(engineered_rows)


if __name__ == "__main__":
    files = glob.glob(os.path.join("data", "*_parsed.csv"))
    for file in files:
        gene = os.path.basename(file).split("_")[0]
        df_parsed = pd.read_csv(file)
        df_engineered = engineer_features(df_parsed)
        out_file = os.path.join("data", f"{gene}_engineered_features.csv")
        df_engineered.to_csv(out_file, index=False)
        print(f"[SUCCESS] {gene}: Engineered physicochemical features for {len(df_engineered)} variants -> {out_file}")
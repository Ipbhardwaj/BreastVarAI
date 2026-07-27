"""
Variant Parser & Data Validation Module for BreastVarAI
Parses HGVS protein notation (e.g., p.Arg175His) into discrete 1-letter amino acid codes
and residue positions, validating raw inputs against malformed records.
"""

import os
import re
import pandas as pd

# Standard 3-letter to 1-letter amino acid dictionary
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
}


def parse_hgvs_protein(title: str):
    """
    Extracts 3-letter amino acids, converts them to 1-letter codes, and retrieves residue position.
    Example: '... (p.Arg175His)' -> ('R', 175, 'H')
    """
    match = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", str(title))
    if match:
        ref_3, pos, alt_3 = match.groups()
        ref_1 = AA_3TO1.get(ref_3)
        alt_1 = AA_3TO1.get(alt_3)
        if ref_1 and alt_1:
            return ref_1, int(pos), alt_1
    return None, None, None


def validate_and_parse_dataset(gene_symbol: str) -> pd.DataFrame:
    """Loads raw ClinVar CSV, parses protein notation, and filters invalid rows."""
    raw_path = os.path.join("data", "raw", f"{gene_symbol}_clinvar.csv")
    if not os.path.exists(raw_path):
        print(f"[WARNING] File not found: {raw_path}")
        return pd.DataFrame()

    df = pd.read_csv(raw_path)

    # Deduplicate ClinVar IDs
    df = df.drop_duplicates(subset=["clinvar_id"]).copy()

    parsed_data = []
    for _, row in df.iterrows():
        ref_aa, pos, alt_aa = parse_hgvs_protein(row["variant_title"])
        if ref_aa and pos and alt_aa:
            parsed_data.append(
                {
                    "clinvar_id": row["clinvar_id"],
                    "gene_symbol": gene_symbol,
                    "variant_title": row["variant_title"],
                    "protein_position": pos,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "clinical_significance": row["clinical_significance"],
                    "review_status": row["review_status"],
                }
            )

    df_parsed = pd.DataFrame(parsed_data)
    out_path = os.path.join("data", f"{gene_symbol}_parsed.csv")
    df_parsed.to_csv(out_path, index=False)
    print(f"[SUCCESS] {gene_symbol}: Validated {len(df_parsed)} / {len(df)} missense variants -> {out_path}")
    return df_parsed


if __name__ == "__main__":
    genes = ["BRCA1", "BRCA2", "TP53", "PALB2", "CHEK2", "ATM"]
    for gene in genes:
        validate_and_parse_dataset(gene)
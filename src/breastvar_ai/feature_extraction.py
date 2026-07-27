"""
Authentic Feature Annotation Module for BreastVarAI
Fetches genuine UniProt functional protein domain features, amino acid positions,
and sequence properties without synthetic or random placeholders.
"""

import os
import re
import pandas as pd
import requests

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# Primary UniProt Accessions for Breast Cancer Target Genes
GENE_TO_UNIPROT = {
    "BRCA1": "P38398",
    "BRCA2": "P51587",
    "TP53": "P04637",
    "PALB2": "Q86YC2",
    "CHEK2": "O96017",
    "ATM": "Q13315",
}


def parse_amino_acid_change(variant_title: str):
    """Extracts position and amino acid substitution from ClinVar variant title (e.g., 'p.Arg175His')."""
    match = re.search(r"p\.([A-Za-z]+)(\d+)([A-Za-z]+)", str(variant_title))
    if match:
        ref_aa, pos, alt_aa = match.groups()
        return ref_aa, int(pos), alt_aa
    return None, None, None


def fetch_uniprot_features(gene_symbol: str) -> dict:
    """Retrieves real protein domain boundaries from the UniProt REST API."""
    uniprot_id = GENE_TO_UNIPROT.get(gene_symbol)
    if not uniprot_id:
        return {}

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        domains = []
        features = data.get("features", [])
        for feat in features:
            if feat.get("type") in ["Domain", "Active site", "Binding site"]:
                location = feat.get("location", {})
                start = location.get("start", {}).get("value")
                end = location.get("end", {}).get("value")
                description = feat.get("description", feat.get("type"))
                if start and end:
                    domains.append({"start": start, "end": end, "name": description})

        return {"uniprot_id": uniprot_id, "domains": domains}
    except Exception as e:
        print(f"[WARNING] Could not fetch UniProt features for {gene_symbol}: {e}")
        return {"uniprot_id": uniprot_id, "domains": []}


def annotate_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Annotates variants with genuine protein positions and UniProt functional domain overlaps."""
    if df.empty:
        return df

    gene_symbol = df["gene_symbol"].iloc[0]
    uniprot_info = fetch_uniprot_features(gene_symbol)
    domains = uniprot_info.get("domains", [])

    annotated_rows = []
    for _, row in df.iterrows():
        ref_aa, pos, alt_aa = parse_amino_acid_change(row["variant_title"])

        # Check if position falls inside a known UniProt functional domain
        in_domain = 0
        domain_name = "None"
        if pos and domains:
            for d in domains:
                if d["start"] <= pos <= d["end"]:
                    in_domain = 1
                    domain_name = d["name"]
                    break

        # Class label: 1 for Pathogenic, 0 for Benign
        target = 1 if "Pathogenic" in row["clinical_significance"] else 0

        annotated_rows.append(
            {
                "clinvar_id": row["clinvar_id"],
                "gene_symbol": gene_symbol,
                "uniprot_id": uniprot_info.get("uniprot_id", ""),
                "variant_title": row["variant_title"],
                "protein_position": pos,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "is_in_functional_domain": in_domain,
                "domain_name": domain_name,
                "clinical_significance": row["clinical_significance"],
                "target_class": target,
            }
        )

    return pd.DataFrame(annotated_rows)


if __name__ == "__main__":
    for gene in GENE_TO_UNIPROT.keys():
        raw_file = os.path.join("data", "raw", f"{gene}_clinvar.csv")
        if os.path.exists(raw_file):
            df_raw = pd.read_csv(raw_file)
            df_annotated = annotate_variants(df_raw)
            processed_path = os.path.join("data", f"{gene}_annotated.csv")
            df_annotated.to_csv(processed_path, index=False)
            print(f"[SUCCESS] Saved authentic UniProt-annotated dataset for {gene} -> {processed_path}")
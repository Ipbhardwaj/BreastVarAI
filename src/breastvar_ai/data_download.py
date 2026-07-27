"""
Authentic Data Collection Module for BreastVarAI
Queries the NCBI ClinVar API via E-utilities to retrieve real missense variants 
for target breast cancer genes.
"""

import os
import time
import pandas as pd
import requests

TARGET_GENES = ["BRCA1", "BRCA2", "TP53", "PALB2", "CHEK2", "ATM"]
RAW_DATA_DIR = os.path.join("data", "raw")
CLINVAR_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def fetch_clinvar_gene_variants(gene_symbol: str, max_results: int = 40) -> pd.DataFrame:
    """Queries NCBI ClinVar E-utilities for authentic missense variants of a target gene."""
    os.makedirs(RAW_DATA_DIR, exist_ok=True)
    print(f"[INFO] Querying NCBI ClinVar API for {gene_symbol} variants...")

    # Working NCBI search term format
    term = f"{gene_symbol}[gene] AND missense variant[variant consequence]"

    search_params = {
        "db": "clinvar",
        "term": term,
        "retmode": "json",
        "retmax": max_results,
    }

    try:
        response = requests.get(CLINVAR_ESEARCH_URL, params=search_params, timeout=15)
        response.raise_for_status()
        search_data = response.json()
        id_list = search_data.get("esearchresult", {}).get("idlist", [])

        if not id_list:
            # Fallback query format
            fallback_term = f"{gene_symbol} missense"
            search_params["term"] = fallback_term
            response = requests.get(CLINVAR_ESEARCH_URL, params=search_params, timeout=15)
            search_data = response.json()
            id_list = search_data.get("esearchresult", {}).get("idlist", [])

        if not id_list:
            print(f"[WARNING] No ClinVar records returned for {gene_symbol}")
            return pd.DataFrame()

        # Fetch summaries for retrieved IDs
        summary_params = {
            "db": "clinvar",
            "id": ",".join(id_list),
            "retmode": "json",
        }

        summary_resp = requests.get(CLINVAR_ESUMMARY_URL, params=summary_params, timeout=15)
        summary_resp.raise_for_status()
        summary_data = summary_resp.json().get("result", {})

        records = []
        for uid in id_list:
            if uid in summary_data:
                item = summary_data[uid]
                title = item.get("title", "")
                germline = item.get("germline_classification", {})
                significance = germline.get("description", "Unknown")

                records.append(
                    {
                        "clinvar_id": uid,
                        "gene_symbol": gene_symbol,
                        "variant_title": title,
                        "clinical_significance": significance,
                        "review_status": germline.get("review_status", "Unreviewed"),
                    }
                )

        df = pd.DataFrame(records)
        output_path = os.path.join(RAW_DATA_DIR, f"{gene_symbol}_clinvar.csv")
        df.to_csv(output_path, index=False)
        print(f"[SUCCESS] Retrieved {len(df)} authentic ClinVar variants for {gene_symbol} -> {output_path}")
        return df

    except Exception as e:
        print(f"[ERROR] Failed to query NCBI ClinVar for {gene_symbol}: {e}")
        return pd.DataFrame()


if __name__ == "__main__":
    for gene in TARGET_GENES:
        fetch_clinvar_gene_variants(gene, max_results=30)
        time.sleep(0.4)
\# BreastVarAI 



> \*\*An Interpretable AI Framework for Predicting the Structural \& Functional Impact of Breast Cancer Variants\*\*



\[!\[Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)

\[!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

\[!\[Tests](https://img.shields.io/badge/pytest-passing-brightgreen.svg)](tests/)



\## Overview

\*\*`BreastVarAI`\*\* is an open-source computational framework designed to predict whether germline and somatic missense variants in hereditary breast cancer genes (\*\*BRCA1, BRCA2, PALB2, TP53, CHEK2, ATM\*\*) disrupt 3D protein structure and function. By integrating 3D structural metrics (AlphaFold/PDB coordinates, solvent accessibility, domain locality) with ensemble machine learning (XGBoost, LightGBM) and TreeSHAP feature attributions, the pipeline classifies variant pathogenicity with complete biophysical interpretability.



\---



\## Workflow Architecture



```text

 ClinVar / COSMIC / gnomAD Variants (BRCA1, BRCA2, PALB2, TP53, CHEK2, ATM)

                                    │

                                    ▼

                3D Structural Mapping (AlphaFold DB / PDB)

                                     │

                                     ▼

          Biophysical Feature Extraction (SASA, Secondary Structure, ΔΔG)

                                   │

                                    ▼

           Ensemble AI Classifiers (XGBoost / LightGBM / Random Forest)

                                    │

                                     ▼

             TreeSHAP Feature Attribution \& Pathogenicity Scores


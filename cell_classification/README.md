## Overview

This directory contains two Jupyter notebooks needed to classify cells. The RF model is developed in the `preprocessing.ipynb` notebook, while the RF model is deployed and prediction probabilities evaluated in the `alignment.ipynb` notebook.

To run the `preprocessing.ipynb` notebook, you need acces to the Excel sheet that defines cell phenotypes based on binary protein expression (`BinaryCellClassification.xlsx`) and the entire raw data set (`WholeTissueConsolidated_data.csv`), which are both available in Zenodo. The main output of this notebook is a .joblib file that contains a trained random forest model that can be used to label cells that do not fit the cell phenotypes defined using protein expression. This file has been precomputed and is also available in Zenodo (`all_nontx_data_rf.joblib`).

To run the `alignment.ipynb` notebook, you need this RF model and two data directories that are zipped and available on Zenodo. These are the metastasis clusters with fluorescence data (`complete_cluster_data.zip`) and the extracted metastases with no fluorescence data (`extracted_met_data.zip`). This notebook first aligns extracted metastases with their associated fluorescence data, then uses the RF model to predict unlabelled cell type labels using this data.

Thus, the structure of this directory should resemble:

```
cell_classification
├── alignment.ipynb
├── all_nontx_data_rf.joblib
├── BinaryCellClassification.xlsx
├── complete_cluster_data/ (from complete_cluster_data.zip)
├── extracted_met_data/ (from extracted_met_data.zip)
├── preprocessing.ipynb
├── README.md
└── WholeTissueConsolidated_data.csv
```

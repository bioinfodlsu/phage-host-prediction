# Phage-Host Interaction Prediction

![badge][badge-jupyter]
![badge][badge-python]
![badge][badge-pandas]
![badge][badge-numpy]
![badge][badge-scipy]
![scikit-learn](https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=flat&logo=scikit-learn&logoColor=white)

[![Actions Status](https://github.com/bioinfodlsu/phage-host-prediction/workflows/Check%20for%20syntax%20errors/badge.svg)](https://github.com/bioinfodlsu/phage-host-prediction/actions)
![badge][badge-github-actions]

## Description
**ABSTRACT**: With the growing interest in using phages to combat antimicrobial resistance, computational methods for predicting phage-host interactions have been explored to help shortlist candidate phages. While existing systems have been successful in integrating multiple features to improve performance, most consider entire proteomes and rely on manual feature engineering, which poses difficulty in selecting the most informative sequence properties to serve as input to the model. In this paper, we focused on the phages' receptor-binding proteins, which are known to be the key machinery for host recognition, and explored different protein language models to automatically encode these protein sequences into meaningful dense embeddings without the need for additional alignment or structural information. Our experiments showed that the use of embeddings of receptor-binding proteins presents improvements over handcrafted genomic and protein sequence features for phage-host interaction prediction. The highest performance was obtained using the transformer-based protein language model ProtT5, resulting in a 3% to 4% increase in the weighted F1 scores across different thresholds for prediction confidence.

## Project Structure
The [`experiments`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments) folder contains the files and scripts for running our model and reproducing our results. Additional (large) requisite files can be downloaded following the instructions provided in the Jupyter notebooks.

### Directories

Directory | Description
-- | --
[`inphared`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/inphared) | Contains the list of phage-host pairs in TSV format. The GenBank and FASTA files with the genomic and protein sequences of the phages, the embeddings of the receptor-binding proteins, and the phage-host-features CSV files should also be saved in this folder
[`preprocessing`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/preprocessing) | Contains text files related to the preprocessing of host information and the selection of annotated RBPs
[`rbp_prediction`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/rbp_prediction) | Contains the JSON file of the trained XGBoost model proposed by [Boeckaerts <i>et al.</i> (2022)](https://www.mdpi.com/1999-4915/14/6/1329) for the computational prediction of RBPs
[`temp`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/temp) | Contains intermediate output files during preprocessng and performance evaluation

### Jupyter Notebooks

### Python Scripts

## Dependencies
The following Python libraries and modules were used. 

Libraries/Modules | Description | License
-- | -- | --
[`regex`](https://github.com/mrabarnett/mrab-regex)
[`nltk`](https://www.nltk.org/)
[`biopython`](https://biopython.org/)
[`ete3`](http://etetoolkit.org/)
[`pandas`](https://pandas.pydata.org/)
[`numpy`](https://numpy.org/)
[`scipy`](https://scipy.org/)
[`scikit-learn`](https://scikit-learn.org/)
[`xgboost`](https://xgboost.readthedocs.io/en/stable/)
[`joblib`](https://joblib.readthedocs.io/en/latest/)
[`numba`](https://numba.pydata.org/)
[`matplotlib`](https://matplotlib.org/)
[`umap-learn`](https://umap-learn.readthedocs.io/en/latest/)

## Authors
- **Mark Edward M. Gonzales** <br>
  mark_gonzales@dlsu.edu.ph 
 
- **Ms. Jennifer C. Ureta** <br>
  jennifer.ureta@dlsu.edu.ph 
  
- **Dr. Anish M.S. Shrestha** <br>
  anish.shrestha@dlsu.edu.ph

This is a research project under the [Bioinformatics Laboratory](https://bioinfodlsu.com/), [Advanced Research Institute for Informatics, Computing and Networking](https://www.dlsu.edu.ph/research/research-centers/adric/), De La Salle University, Philippines.

[badge-jupyter]: https://img.shields.io/badge/Jupyter-F37626.svg?&style=flat&logo=Jupyter&logoColor=white
[badge-python]: https://img.shields.io/badge/python-3670A0?style=flat&logo=python&logoColor=white
[badge-pandas]: https://img.shields.io/badge/Pandas-2C2D72?style=flat&logo=pandas&logoColor=white
[badge-numpy]: https://img.shields.io/badge/Numpy-777BB4?style=flat&logo=numpy&logoColor=white
[badge-scipy]: https://img.shields.io/badge/SciPy-654FF0?style=flat&logo=SciPy&logoColor=white
[badge-github-actions]: https://img.shields.io/badge/GitHub_Actions-2088FF?style=flat&logo=github-actions&logoColor=white

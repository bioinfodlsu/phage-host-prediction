# Phage-Host Interaction Prediction

![badge][badge-jupyter]
![badge][badge-python]
![badge][badge-pandas]
![badge][badge-numpy]
![badge][badge-scipy]
![scikit-learn](https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=flat&logo=scikit-learn&logoColor=white)  <br>
[![Actions Status](https://github.com/bioinfodlsu/phage-host-prediction/workflows/Check%20for%20syntax%20errors/badge.svg)](https://github.com/bioinfodlsu/phage-host-prediction/actions)
![badge][badge-github-actions]

Our paper can be accessed via this [link]().

## Table of Contents
- [Description](https://github.com/bioinfodlsu/phage-host-prediction#description)
- [Project Structure](https://github.com/bioinfodlsu/phage-host-prediction#project-structure)
  - [Directories](https://github.com/bioinfodlsu/phage-host-prediction#directories)
  - [Jupyter Notebooks](https://github.com/bioinfodlsu/phage-host-prediction#jupyter-notebooks)
  - [Python Scripts](https://github.com/bioinfodlsu/phage-host-prediction#python-scripts)
  - [Folder Structure](https://github.com/bioinfodlsu/phage-host-prediction#folder-structure)
- [Dependencies](https://github.com/bioinfodlsu/phage-host-prediction#dependencies)
- [Authors](https://github.com/bioinfodlsu/phage-host-prediction#authors)

## Description
**ABSTRACT**: With the growing interest in using phages to combat antimicrobial resistance, computational methods for predicting phage-host interactions have been explored to help shortlist candidate phages. While existing systems have been successful in integrating multiple features to improve performance, most consider entire proteomes and rely on manual feature engineering, which poses difficulty in selecting the most informative sequence properties to serve as input to the model. In this paper, we focused on the phages' receptor-binding proteins, which are known to be the key machinery for host recognition, and explored different protein language models to automatically encode these protein sequences into meaningful dense embeddings without the need for additional alignment or structural information. Our experiments showed that the use of embeddings of receptor-binding proteins presents improvements over handcrafted genomic and protein sequence features for phage-host interaction prediction. The highest performance was obtained using the transformer-based protein language model ProtT5, resulting in a 3% to 4% increase in the weighted F1 scores across different thresholds for prediction confidence.

<img src="https://github.com/bioinfodlsu/phage-host-prediction/blob/main/figure.png?raw=True" alt="Teaser Figure" width = 800> 

↑ *Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction#table-of-contents).*

## Project Structure
The [`experiments`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments) folder contains the files and scripts for running our model and reproducing our results. Note that additional (large) files have to be downloaded (or generated) following the instructions in the Jupyter notebooks.

### Directories

Directory | Description
-- | --
[`inphared`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/inphared) | Contains the list of phage-host pairs in TSV format. The GenBank and FASTA files with the genomic and protein sequences of the phages, the embeddings of the receptor-binding proteins, and the phage-host-features CSV files should also be saved in this folder
[`preprocessing`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/preprocessing) | Contains text files related to the preprocessing of host information and the selection of annotated RBPs
[`rbp_prediction`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/rbp_prediction) | Contains the JSON file of the trained XGBoost model proposed by [Boeckaerts <i>et al.</i> (2022)](https://www.mdpi.com/1999-4915/14/6/1329) for the computational prediction of receptor-binding proteins. Downloaded from this [repository](https://github.com/dimiboeckaerts/PhageRBPdetection/blob/main/data/RBPdetect_xgb_model.json) (under the MIT License)
[`temp`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/temp) | Contains intermediate output files during preprocessng and performance evaluation

↑ *Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction#table-of-contents).*

### Jupyter Notebooks
Each notebook provides detailed instructions related to the required and output files, including the download links and where to save them.

Notebook | Description | Required Files | Output Files
-- | -- | -- | --
[`1. Sequence Preprocessing.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/1.%20Sequence%20Preprocessing.ipynb) | Preprocessing of host information and selection of annotated receptor-binding proteins | GenomesDB ([Partial](https://millardlab-inphared.s3.climb.ac.uk/GenomesDB_20201412.tar.gz). Complete populating following the instructions in the notebook), <br> [GenBank file of phage genomes and/or proteomes](https://drive.google.com/file/d/14LG1iGa1CqPbAjofZT1EY8VKnE8Iy45Q/view?usp=sharing) | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing)
[`2. Exploratory Data Analysis.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/2.%20Exploratory%20Data%20Analysis.ipynb) | Exploratory data analysis | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)), <br> [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing) | &ndash;
[`3. RBP Computational Prediction.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/3.%20RBP%20Computational%20Prediction.ipynb) | Computational prediction of receptor-binding proteins | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))
[`4. Protein Embedding Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/4.%20Protein%20Embedding%20Generation.ipynb) | Generation of protein embeddings | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing) | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))
[`5. Data Consolidation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/5.%20Data%20Consolidation.ipynb) | Generation of phage-host-features CSV files | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing), <br> Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)
[`6. Classifier Building & Evaluation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/6.%20Classifier%20Building%20%26%20Evaluation.ipynb) | Construction of phage-host interaction model and performance evaluation | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing) | [Trained models](https://drive.google.com/drive/folders/1U5ugmkhD4LHElYnLj3B8Xt2TcPx-TOjB?usp=sharing)
[`7. Visualization.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/7.%20Visualization.ipynb) | Plotting of <i>t</i>-SNE and UMAP projections | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing) | &ndash;

↑ *Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction#table-of-contents).*

### Python Scripts

Script | Description |
-- | --
[`ClassificationUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/ClassificationUtil.py) | Contains the utility functions for the generation of the phage-host-features CSV files, construction of the phage-host interaction model, and performance evaluation
[`ConstantsUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/ConstantsUtil.py) | Contains the constants used in the notebooks and scripts
[`EDAUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/EDAUtil.py) | Contains the utility functions for exploratory data analysis
[`RBPPredictionUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/RBPPredictionUtil.py) | Contains the utility functions for the computational prediction of receptor-binding proteins
[`SequenceParsing.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/SequenceParsingUtil.py) | Contains the utility functions for preprocessing host information and selecting annotated receptor-binding proteins
[`boeckaerts.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/boeckaerts.py) | Contains the utility functions written by [Boeckaerts <i>et al.</i> (2021)](https://www.nature.com/articles/s41598-021-81063-4) for running his phage-host interaction prediction tool (with which we benchmarked our model). Downloaded from this [repository](https://github.com/dimiboeckaerts/BacteriophageHostPrediction/blob/master/RBP_functions.py) (under the MIT License)

↑ *Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction#table-of-contents).*

### Folder Structure
Once you have cloned this repository and finished downloading (or generating) all the additional required files following the instructions in the Jupyter notebooks, your folder structure should be similar to the one below:

- `phage-host-prediction` (root)
  - `datasets` 
    - `inphared`
      - `inphared`
        - `GenomesDB` (Downoad [partial](https://millardlab-inphared.s3.climb.ac.uk/GenomesDB_20201412.tar.gz). Complete populating following the instructions [here](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/1.%20Sequence%20Preprocessing.ipynb))
          - `AB002632`
          - ...
  - `experiments`
    - `inphared`
      - `data` ([Download](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing))
        - `rbp.csv`
        - `rbp_embeddings_esm.csv`
        - ...
      - `embeddings` (Download [Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))
        - `esm`
        - `esm1b`
        - ...
      - `fasta` ([Download](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing))
        - `hypothetical`
        - `nucleotide`
        - `rbp`
      - `16Sep2022_data_excluding_refseq.tsv`
      - `16Sep2022_phages_downloaded_from_genbank.gb` ([Download](https://drive.google.com/file/d/14LG1iGa1CqPbAjofZT1EY8VKnE8Iy45Q/view?usp=sharing))
    - `models` ([Download](https://drive.google.com/drive/folders/1U5ugmkhD4LHElYnLj3B8Xt2TcPx-TOjB?usp=sharing))
      - `boeckaerts.joblib`
      - `esm.joblib`
      - ...
    - `preprocessing`
    - `rbp_prediction`
    - `temp`
    - `1. Sequence Preprocessing.ipynb`
    - ...
    - `ClassificationUtil.py`
    - ...

↑ *Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction#table-of-contents).*

## Dependencies
The following Python libraries and modules were used: 

Libraries/Modules | Description | License
-- | -- | --
[`regex`](https://github.com/mrabarnett/mrab-regex) | Provides additional functionality over the standard [`re`](https://docs.python.org/3/library/re.html) module while maintaining backwards-compatibility	| Apache License 2.0 
[`nltk`](https://www.nltk.org/) | Provides interfaces to corpora and lexical resources, along with a suite of text processing libraries for classification, tokenization, stemming, tagging, parsing, and semantic reasoning	| Apache License 2.0
[`biopython`](https://biopython.org/) | Provides tools for computational molecular biology | Biopython License Agreement, BSD 3-Clause License
[`ete3`](http://etetoolkit.org/) | Provides functions for automated manipulation, analysis, and visualization of phylogenetic trees | GNU General Public License v3.0
[`pandas`](https://pandas.pydata.org/) | Provides functions for data analysis and manipulation	| BSD 3-Clause "New" or "Revised" License
[`numpy`](https://numpy.org/) | Provides a multidimensional array object, various derived objects, and an assortment of routines for fast operations on arrays | BSD 3-Clause "New" or "Revised" License 
[`scipy`](https://scipy.org/) | Provides efficient numerical routines, such as those for numerical integration, interpolation, optimization, linear algebra, and statistics | BSD 3-Clause "New" or "Revised" License
[`scikit-learn`](https://scikit-learn.org/) | Provides efficient tools for predictive data analysis | BSD 3-Clause "New" or "Revised" License
[`xgboost`](https://xgboost.readthedocs.io/en/stable/) | Implements machine learning algorithms under the gradient boosting framework | Apache License 2.0 
[`joblib`](https://joblib.readthedocs.io/en/latest/) | Provides tools for lightweight pipelining in Python | BSD 3-Clause "New" or "Revised" License
[`numba`](https://numba.pydata.org/) | Translates Python functions to optimized machine code at runtime using the industry-standard LLVM compiler library | BSD 2-Clause "Simplified" License
[`matplotlib`](https://matplotlib.org/) | Provides functions for creating static, animated, and interactive visualizations | Matplotlib License (BSD-Compatible)
[`umap-learn`](https://umap-learn.readthedocs.io/en/latest/) | Implements uniform manifold approximation and projection, a dimension reduction technique that can be used for visualisation similarly and general non-linear dimension reduction | BSD 3-Clause "New" or "Revised" License

*The descriptions are taken from their respective websites.*

↑ *Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction#table-of-contents).*

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

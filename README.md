# Protein embeddings improve phage-host interaction prediction

![badge][badge-jupyter]
![badge][badge-python]
![badge][badge-pandas]
![badge][badge-numpy]
![badge][badge-scipy]
![scikit-learn](https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=flat&logo=scikit-learn&logoColor=white)  <br>
[![Actions Status](https://github.com/bioinfodlsu/phage-host-prediction/workflows/Check%20for%20syntax%20errors/badge.svg)](https://github.com/bioinfodlsu/phage-host-prediction/actions)
![badge][badge-github-actions]

**This work was accepted for publication in *PLOS ONE.***
- The final version of our paper (as published in *PLOS ONE*) can be accessed via this [link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0289030).
- Presenting this paper, the lead author (Mark Edward M. Gonzales) won 2nd Prize at the [2023 Magsaysay Future Engineers/Technologists Award](https://fb.watch/oKx7G6gwLi/). This award is conferred by the National Academy of Science and Technology, the highest recognition and scientific advisory body of the Philippines, to recognize outstanding research outputs on engineering and technology at the collegiate level.

If you find our work useful, please consider citing:

```
@article{10.1371/journal.pone.0289030,
    doi = {10.1371/journal.pone.0289030},
    author = {Gonzales, Mark Edward M. AND Ureta, Jennifer C. AND Shrestha, Anish M. S.},
    journal = {PLOS ONE},
    publisher = {Public Library of Science},
    title = {Protein embeddings improve phage-host interaction prediction},
    year = {2023},
    month = {07},
    volume = {18},
    url = {https://doi.org/10.1371/journal.pone.0289030},
    pages = {1-22},
    number = {7}
}
```

## Table of Contents
- [Description](https://github.com/bioinfodlsu/phage-host-prediction#description)
- [Project Structure](https://github.com/bioinfodlsu/phage-host-prediction#project-structure)
  - [Directories](https://github.com/bioinfodlsu/phage-host-prediction#directories)
  - [Jupyter Notebooks](https://github.com/bioinfodlsu/phage-host-prediction#jupyter-notebooks)
  - [Python Scripts](https://github.com/bioinfodlsu/phage-host-prediction#python-scripts)
  - [Folder Structure](https://github.com/bioinfodlsu/phage-host-prediction#folder-structure)
- [Environment & Dependencies](https://github.com/bioinfodlsu/phage-host-prediction#environment--dependencies)
- [Authors](https://github.com/bioinfodlsu/phage-host-prediction#authors)

## Description
**ABSTRACT**: With the growing interest in using phages to combat antimicrobial resistance, computational methods for predicting phage-host interactions have been explored to help shortlist candidate phages. Most existing models consider entire proteomes and rely on manual feature engineering, which poses difficulty in selecting the most informative sequence properties to serve as input to the model. In this paper, we framed phage-host interaction prediction as a multiclass classification problem that takes as input the embeddings of a phage's receptor-binding proteins, which are known to be the key machinery for host recognition, and predicts the host genus. We explored different protein language models to automatically encode these protein sequences into dense embeddings without the need for additional alignment or structural information. We show that the use of embeddings of receptor-binding proteins presents improvements over handcrafted genomic and protein sequence features. The highest performance was obtained using the transformer-based protein language model ProtT5, resulting in a 3% to 4% increase in weighted F1 and recall scores across different prediction confidence thresholds, compared to using selected handcrafted sequence features.

<img src="https://github.com/bioinfodlsu/phage-host-prediction/blob/main/figure.png?raw=True" alt="Teaser Figure" width = 800> 

**AUTHOR SUMMARY**: Antimicrobial resistance is among the major global health issues at present. As alternatives to the usual antibiotics, drug formulations based on phages (bacteria-infecting viruses) have received increased interest, as phages are known to attack only a narrow range of bacterial hosts and antagonize the target pathogen with minimal side effects. The screening of candidate phages has recently been facilitated through the use of machine learning models for inferring phage-host pairs. The performance of these models relies heavily on the transformation of raw biological sequences into a collection of numerical features. However, since a wide array of potentially informative features can be extracted from sequences, selecting the most relevant ones is challenging. Our approach eliminates the need for this manual feature engineering by employing protein language models to automatically generate numerical representations for specific subsets of tail proteins known as receptor-binding proteins. These proteins are responsible for a phage's initial contact with the host bacterium and are thus regarded as important determinants of host specificity. Our results show that this approach presents improvements over using handcrafted genomic and protein sequence features in predicting phage-host interaction.

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
[`3.1. RBP FASTA Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/3.1.%20RBP%20FASTA%20Generation.ipynb) | Generation of the FASTA files containing the RBP protein sequences] | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing)
[`4. Protein Embedding Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/4.%20Protein%20Embedding%20Generation.ipynb) | Generation of protein embeddings | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing) | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))
[`5. Data Consolidation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/5.%20Data%20Consolidation.ipynb) | Generation of phage-host-features CSV files | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing), <br> Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)
[`6. Classifier Building & Evaluation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/6.%20Classifier%20Building%20%26%20Evaluation.ipynb) | Construction of phage-host interaction model and performance evaluation | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing) | [Trained models](https://drive.google.com/drive/folders/1U5ugmkhD4LHElYnLj3B8Xt2TcPx-TOjB?usp=sharing)
[`6.1. Additional Model Evaluation (Specificity + PR Curve).ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/6.1.%20Additional%20Model%20Evaluation%20(Specificity%20%2B%20PR%20Curve).ipynb) | Addition of metrics for model evaluation | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)| &ndash;
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
[`boeckaerts.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/boeckaerts.py) | Contains the utility functions written by [Boeckaerts <i>et al.</i> (2021)](https://www.nature.com/articles/s41598-021-81063-4) for running their phage-host interaction prediction tool (with which we benchmarked our model). Downloaded from this [repository](https://github.com/dimiboeckaerts/BacteriophageHostPrediction/blob/master/RBP_functions.py) (under the MIT License)

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
        - `complete` 
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

## Environment & Dependencies

⚠️ **UPDATE (06/12/2023):** In May 2023, Google Colab switched its default runtime to Python 3.10. However, one of our project's dependencies, [`bio-embeddings`](https://docs.bioembeddings.com/v0.2.3/) (v0.2.3), seems to be incompatible with Python 3.10. 

If the memory requirement of loading pretrained protein language models ([`4. Protein Embedding Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/4.%20Protein%20Embedding%20Generation.ipynb)) is too heavy for your local machine, an alternative cloud-based service with GPU is Paperspace; you may try using either its [PyTorch 1.12 runtime](https://docs.paperspace.com/gradient/notebooks/runtimes/#recommended-runtimes) (which, as of writing, uses Python 3.9) or [Python 3.9 runtime](https://docs.paperspace.com/gradient/notebooks/runtimes/#previous-runtime-versions).

### Operating System

One of our project's dependencies, [`bio_embeddings`](https://docs.bioembeddings.com/v0.2.3/), was developed for Unix and Unix-like operating systems. If you are running this project on Windows, consider using [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) (WSL) or a virtual machine.

### Dependencies

We recommend using Python 3.9 to run this project. Thanks to Dr. Paul K. Yu (paul_yu_a@dlsu.edu.ph) for sharing his environment configuration ([`environment.yaml`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/environment.yaml)).

The dependencies can be installed via [Conda](https://docs.conda.io/en/latest/), an open-source package and environment management system. Run the following command to create a virtual environment with the dependencies installed:
```
conda env create -f environment.yaml
```

To activate this environment, run the following command:
```
conda activate phage-host-prediction
```

The complete list of Python libraries and modules used in this project (excluding those that are part of the Python Standard Library) is as follows: 

Library/Module | Description | License
-- | -- | --
[`regex`](https://github.com/mrabarnett/mrab-regex) | Provides additional functionality over the standard [`re`](https://docs.python.org/3/library/re.html) module while maintaining backwards-compatibility	| Apache License 2.0 
[`nltk`](https://www.nltk.org/) | Provides interfaces to corpora and lexical resources, along with a suite of text processing libraries for classification, tokenization, stemming, tagging, parsing, and semantic reasoning	| Apache License 2.0
[`biopython`](https://biopython.org/) | Provides tools for computational molecular biology | Biopython License Agreement, BSD 3-Clause License
[`ete3`](http://etetoolkit.org/) | Provides functions for automated manipulation, analysis, and visualization of phylogenetic trees | GNU General Public License v3.0
[`pandas`](https://pandas.pydata.org/) | Provides functions for data analysis and manipulation	| BSD 3-Clause "New" or "Revised" License
[`numpy`](https://numpy.org/) | Provides a multidimensional array object, various derived objects, and an assortment of routines for fast operations on arrays | BSD 3-Clause "New" or "Revised" License 
[`scipy`](https://scipy.org/) | Provides efficient numerical routines, such as those for numerical integration, interpolation, optimization, linear algebra, and statistics | BSD 3-Clause "New" or "Revised" License
[`scikit-learn`](https://scikit-learn.org/) | Provides efficient tools for predictive data analysis | BSD 3-Clause "New" or "Revised" License
[`imbalanced-learn`](https://imbalanced-learn.org/stable/) | Provides tools when dealing with classification with imbalanced classes | MIT License
[`pyyaml`](https://pyyaml.org/) | Supports standard YAML tags and provides Python-specific tags that allow to represent an arbitrary Python object | MIT License
[`xgboost`](https://xgboost.readthedocs.io/en/stable/) | Implements machine learning algorithms under the gradient boosting framework | Apache License 2.0 
[`joblib`](https://joblib.readthedocs.io/en/latest/) | Provides tools for lightweight pipelining in Python | BSD 3-Clause "New" or "Revised" License
[`numba`](https://numba.pydata.org/) | Translates Python functions to optimized machine code at runtime using the industry-standard LLVM compiler library | BSD 2-Clause "Simplified" License
[`matplotlib`](https://matplotlib.org/) | Provides functions for creating static, animated, and interactive visualizations | Matplotlib License (BSD-Compatible)
[`jsonnet`](https://jsonnet.org/) | Domain-specific language for JSON | Apache License 2.0 
[`cudatoolkit`](https://developer.nvidia.com/cuda-toolkit) | Parallel computing platform and programming model for general computing on GPUs | NVIDIA Software License
[`bio-embeddings`](https://docs.bioembeddings.com/v0.2.3/) | Provides an interface for the use of language model-based biological sequence representations for transfer-learning | MIT License
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

This research was partly funded by the [Department of Science and Technology &ndash; Philippine Council for Health Research and Development](https://www.pchrd.dost.gov.ph/) (DOST-PCHRD) under the [e-Asia JRP 2021 Alternative therapeutics to tackle AMR pathogens (ATTACK-AMR) program](https://www.the-easia.org/jrp/projects/project_76.html). The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

[badge-jupyter]: https://img.shields.io/badge/Jupyter-F37626.svg?&style=flat&logo=Jupyter&logoColor=white
[badge-python]: https://img.shields.io/badge/python-3670A0?style=flat&logo=python&logoColor=white
[badge-pandas]: https://img.shields.io/badge/Pandas-2C2D72?style=flat&logo=pandas&logoColor=white
[badge-numpy]: https://img.shields.io/badge/Numpy-777BB4?style=flat&logo=numpy&logoColor=white
[badge-scipy]: https://img.shields.io/badge/SciPy-654FF0?style=flat&logo=SciPy&logoColor=white
[badge-github-actions]: https://img.shields.io/badge/GitHub_Actions-2088FF?style=flat&logo=github-actions&logoColor=white

# PHIEmbed (Phage-Host Interaction Prediction with Protein Embeddings)

![badge][badge-jupyter]
![badge][badge-python]
![badge][badge-pandas]
![badge][badge-numpy]
![badge][badge-scipy]
![scikit-learn](https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=flat&logo=scikit-learn&logoColor=white) <br>
[![Actions Status](https://github.com/bioinfodlsu/phage-host-prediction/workflows/Check%20for%20syntax%20errors/badge.svg)](https://github.com/bioinfodlsu/phage-host-prediction/actions)
[![Actions Status](https://github.com/bioinfodlsu/phage-host-prediction/workflows/Run%20Black%20formatter/badge.svg)](https://github.com/bioinfodlsu/phage-host-prediction/actions)
![badge][badge-github-actions]

**PHIEmbed** is a phage-host interaction prediction tool that uses protein language models to represent the receptor-binding proteins of phages. It presents improvements over using handcrafted (manually feature-engineered) sequence properties and eliminates the need to manually extract and select features from phage sequences.

**Paper**: https://doi.org/10.1371/journal.pone.0289030

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

You can also find PHIEmbed on [bio.tools](https://bio.tools/phiembed).

## Table of Contents

-   [📰 News](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#-news)
-   [♾️ Run on Google Colab](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#%EF%B8%8F-run-on-google-colab)
-   [🚀 Installation & Usage](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#-installation--usage)
-   [📚 Description](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#-description)
-   [🧪 Reproducing Our Results](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#-reproducing-our-results)
-   [💻 Authors](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#-authors)

## 📰 News

-   **1 Sep 2024** - We created a Google Colab notebook for running our tool. Instructions [here](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#%EF%B8%8F-run-on-google-colab).
  
-   **24 Apr 2024** - We added scripts to simplify running and training our tool. Instructions [here](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#-installation--usage).

-   **23 Feb 2024** - We presented our work at the **eAsia AMR Workshop 2024** held virtually and in person in Tokyo, Japan, and attended by antimicrobial resistance (AMR) researchers from Thailand, USA, Australia, Japan, and the Philippines. Slides [here](https://docs.google.com/presentation/d/1rnMAg5fIVFuK5JxIQQOAh5311GYgiOOeVVvdttcRY6I/edit?usp=sharing).

-   **01 Dec 2023** - Presenting this work, the lead author (Mark Edward M. Gonzales) won **2nd Prize at the 2023 Magsaysay Future Engineers/Technologists Award**. This award is conferred by the National Academy of Science and Technology, the highest recognition and scientific advisory body of the Philippines, to recognize outstanding research outputs on engineering and technology at the collegiate level. Presentation [here](https://fb.watch/oKx7G6gwLi/) (29:35&ndash;39:51) and slides [here](https://docs.google.com/presentation/d/1Rdjy6l3gnIzcRnAccq2sddltEdJV4g-V1frYc4ge4tQ/edit?usp=sharing).

-   **07 Jul 2023** - Our [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0289030) was accepted for publication in _**PLOS ONE**_.

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

## ♾️ Run on Google Colab

[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://phiembed.bioinfodlsu.com/)

**You can readily run PHIEmbed on Google Colab, without the need to install anything on your own computer: [http://phiembed.bioinfodlsu.com](http://phiembed.bioinfodlsu.com)**

## 🚀 Installation & Usage

**Operating System**: Windows, Linux, or macOS <br>

Clone the repository:

```
git clone https://github.com/bioinfodlsu/phage-host-prediction
cd phage-host-prediction
```

Create a virtual environment with all the necessary dependencies installed via Conda (we recommend using [Miniconda](https://docs.anaconda.com/free/miniconda/index.html)):

```
conda env create -f environment.yaml
```

Activate this environment by running:

```
conda activate PHIEmbed
```

### Running PHIEmbed

```
python3 phiembed.py --input <input_fasta> --model <model_joblib> --output <results_dir>
```

-   Replace `<input_fasta>` with the path to the FASTA file containing the receptor-binding protein sequences. A sample FASTA file is provided [here](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/sample.fasta).
-   Replace `<model_joblib>` with the path to the trained model (recognized format: joblib or compressed joblib, framework: scikit-learn). Download our trained model from this [link](https://drive.google.com/file/d/1bRloKMtPnp8QTOHx5IvSx_-8BspdVKNQ/view?usp=sharing). No need to uncompress, but doing so will speed up loading the model albeit at the cost of additional storage requirements. Refer to this [guide](https://joblib.readthedocs.io/en/latest/generated/joblib.dump.html) for the list of accepted compressed formats.
-   Replace `<results_dir>` with the path to the directory to which the results of running PHIEmbed will be written. The results of running PHIEmbed on the sample FASTA file are provided [here](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/sample_results).

The results for each protein are written to a CSV file (without a header row). Each row contains two comma-separated values: a host genus and the corresponding prediction score (class probability). The rows are sorted in order of decreasing prediction score. Hence, the first row pertains to the top-ranked prediction.

Under the hood, this script first converts each sequence into a protein embedding using ProtT5 (the top-performing protein language model based on our experiments) and then passes the embedding to a random forest classifier trained on our entire [dataset](https://drive.google.com/file/d/1icEenU5Sv-7i9pUycaQfNC1Imhrg3sEN/view?usp=sharing). If your machine has a GPU, it will automatically be used to accelerate the protein embedding generation step.

**Note:** Running this script for the first time may take a few extra minutes since it involves downloading a model (ProtT5, around 2 GB) from Hugging Face.

### Training PHIEmbed

```
python3 train.py --input <training_dataset>
```

-   Replace `<training_dataset>` with the path to the training dataset. A sample can be downloaded [here](https://drive.google.com/file/d/1icEenU5Sv-7i9pUycaQfNC1Imhrg3sEN/view?usp=sharing).
-   The number of threads to be used for training can be specified using `--threads`. By default, it is set to -1 (that is, all threads are to be used).

The training dataset should be formatted as a CSV file (without a header row) where each row corresponds to a training sample. The first column is for the protein IDs, the second column is for the host genera, and the next 1,024 columns are for the components of the ProtT5 embeddings.

This script will output a gzip-compressed, serialized version of the trained model with filename `phiembed_trained.joblib.gz`.

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

## 📚 Description

**Motivation**: With the growing interest in using phages to combat antimicrobial resistance, computational methods for predicting phage-host interactions have been explored to help shortlist candidate phages. Most existing models consider entire proteomes and rely on manual feature engineering, which poses difficulty in selecting the most informative sequence properties to serve as input to the model.

**Method**: In this paper, we framed phage-host interaction prediction as a multiclass classification problem that takes as input the embeddings of a phage's receptor-binding proteins, which are known to be the key machinery for host recognition, and predicts the host genus. We explored different protein language models to automatically encode these protein sequences into dense embeddings without the need for additional alignment or structural information.

**Results**: We show that the use of embeddings of receptor-binding proteins presents improvements over handcrafted genomic and protein sequence features. The highest performance was obtained using the transformer-based protein language model ProtT5, resulting in a 3% to 4% increase in weighted F1 and recall scores across different prediction confidence thresholds, compared to using selected handcrafted sequence features.

<img src="https://github.com/bioinfodlsu/phage-host-prediction/blob/main/figure.png?raw=True" alt="Teaser Figure" width = 800>

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

## 🧪 Reproducing Our Results

### Project Structure

The [`experiments`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments) folder contains the files and scripts for reproducing our results. Note that additional (large) files have to be downloaded (or generated) following the instructions in the Jupyter notebooks.

<details>
  <summary>Click here to show/hide the list of directories, Jupyter notebooks, and Python scripts, as well as the folder structure.</summary>

#### Directories

| Directory                                                                                                     | Description                                                                                                                                                                                                                                                                                                                                                       |
| ------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`inphared`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/inphared)             | Contains the list of phage-host pairs in TSV format. The GenBank and FASTA files with the genomic and protein sequences of the phages, the embeddings of the receptor-binding proteins, and the phage-host-features CSV files should also be saved in this folder                                                                                                 |
| [`preprocessing`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/preprocessing)   | Contains text files related to the preprocessing of host information and the selection of annotated RBPs                                                                                                                                                                                                                                                          |
| [`rbp_prediction`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/rbp_prediction) | Contains the JSON file of the trained XGBoost model proposed by [Boeckaerts <i>et al.</i> (2022)](https://www.mdpi.com/1999-4915/14/6/1329) for the computational prediction of receptor-binding proteins. Downloaded from this [repository](https://github.com/dimiboeckaerts/PhageRBPdetection/blob/main/data/RBPdetect_xgb_model.json) (under the MIT License) |
| [`temp`](https://github.com/bioinfodlsu/phage-host-prediction/tree/main/experiments/temp)                     | Contains intermediate output files during preprocessing and performance evaluation                                                                                                                                                                                                                                                                                |

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

#### Jupyter Notebooks

Each notebook provides detailed instructions related to the required and output files, including the download links and where to save them.

| Notebook                                                                                                                                                                                                                          | Description                                                                            | Required Files                                                                                                                                                                                                                                                                                                                                                    | Output Files                                                                                                                                                                                                           |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`1. Sequence Preprocessing.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/1.%20Sequence%20Preprocessing.ipynb)                                                                               | Preprocessing of host information and selection of annotated receptor-binding proteins | GenomesDB ([Partial](https://millardlab-inphared.s3.climb.ac.uk/GenomesDB_20201412.tar.gz). Complete populating following the instructions in the notebook), <br> [GenBank file of phage genomes and/or proteomes](https://drive.google.com/file/d/14LG1iGa1CqPbAjofZT1EY8VKnE8Iy45Q/view?usp=sharing)                                                            | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing)                                                                                   |
| [`2. Exploratory Data Analysis.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/2.%20Exploratory%20Data%20Analysis.ipynb)                                                                       | Exploratory data analysis                                                              | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)), <br> [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)                | &ndash;                                                                                                                                                                                                                |
| [`3. RBP Computational Prediction.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/3.%20RBP%20Computational%20Prediction.ipynb)                                                                 | Computational prediction of receptor-binding proteins                                  | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))                                                                                                                                            | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) |
| [`3.1. RBP FASTA Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/3.1.%20RBP%20FASTA%20Generation.ipynb)                                                                             | Generation of the FASTA files containing the RBP protein sequences                     | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))                                                                                                                                            | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing)                                                                                   |
| [`4. Protein Embedding Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/4.%20Protein%20Embedding%20Generation.ipynb)                                                                 | Generation of protein embeddings                                                       | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing)                                                                                                                                                                                                                              | Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) |
| [`5. Data Consolidation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/5.%20Data%20Consolidation.ipynb)                                                                                       | Generation of phage-host-features CSV files                                            | [FASTA files of genomic and protein sequences](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing), <br> Protein embeddings ([Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing)) | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)                                                                                                  |
| [`6. Classifier Building & Evaluation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/6.%20Classifier%20Building%20%26%20Evaluation.ipynb)                                                     | Construction of phage-host interaction model and performance evaluation                | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)                                                                                                                                                                                                                                             | [Trained models](https://drive.google.com/drive/folders/1U5ugmkhD4LHElYnLj3B8Xt2TcPx-TOjB?usp=sharing)                                                                                                                 |
| [`6.1. Additional Model Evaluation (Specificity + PR Curve).ipynb`](<https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/6.1.%20Additional%20Model%20Evaluation%20(Specificity%20%2B%20PR%20Curve).ipynb>) | Addition of metrics for model evaluation                                               | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)                                                                                                                                                                                                                                             | &ndash;                                                                                                                                                                                                                |
| [`7. Visualization.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/7.%20Visualization.ipynb)                                                                                                   | Plotting of <i>t</i>-SNE and UMAP projections                                          | [Phage-host-features CSV files](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing)                                                                                                                                                                                                                                             | &ndash;                                                                                                                                                                                                                |

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

#### Python Scripts

| Script                                                                                                                        | Description                                                                                                                                                                                                                                                                                                                                                                           |
| ----------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`ClassificationUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/ClassificationUtil.py)   | Contains the utility functions for the generation of the phage-host-features CSV files, construction of the phage-host interaction model, and performance evaluation                                                                                                                                                                                                                  |
| [`ConstantsUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/ConstantsUtil.py)             | Contains the constants used in the notebooks and scripts                                                                                                                                                                                                                                                                                                                              |
| [`EDAUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/EDAUtil.py)                         | Contains the utility functions for exploratory data analysis                                                                                                                                                                                                                                                                                                                          |
| [`RBPPredictionUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/RBPPredictionUtil.py)     | Contains the utility functions for the computational prediction of receptor-binding proteins                                                                                                                                                                                                                                                                                          |
| [`SequenceParsingUtil.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/SequenceParsingUtil.py) | Contains the utility functions for preprocessing host information and selecting annotated receptor-binding proteins                                                                                                                                                                                                                                                                   |
| [`boeckaerts.py`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/boeckaerts.py)                   | Contains the utility functions written by [Boeckaerts <i>et al.</i> (2021)](https://www.nature.com/articles/s41598-021-81063-4) for running their phage-host interaction prediction tool (with which we benchmarked our model). Downloaded from this [repository](https://github.com/dimiboeckaerts/BacteriophageHostPrediction/blob/master/RBP_functions.py) (under the MIT License) |

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

#### Folder Structure

Once you have cloned this repository and finished downloading (or generating) all the additional required files following the instructions in the Jupyter notebooks, your folder structure should be similar to the one below:

-   `phage-host-prediction` (root)
    -   `datasets`
        -   `inphared`
            -   `inphared`
                -   `GenomesDB` (Downoad [partial](https://millardlab-inphared.s3.climb.ac.uk/GenomesDB_20201412.tar.gz). Complete populating following the instructions [here](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/1.%20Sequence%20Preprocessing.ipynb))
                    -   `AB002632`
                    -   ...
    -   `experiments`
        -   `inphared`
            -   `data` ([Download](https://drive.google.com/drive/folders/1xNoA6dxkN4jzVNCg_7YNjdPZzl51Jo9M?usp=sharing))
                -   `rbp.csv`
                -   `rbp_embeddings_esm.csv`
                -   ...
            -   `embeddings` (Download [Part 1](https://drive.google.com/drive/folders/1deenrDQIr3xcl9QCYH-nPhmpY8x2drQw?usp=sharing) and [Part 2](https://drive.google.com/drive/folders/1jnBFNsC6zJISkc6IAz56257MSXKjY0Ez?usp=sharing))
                -   `esm`
                -   `esm1b`
                -   ...
            -   `fasta` ([Download](https://drive.google.com/drive/folders/16ZBXZCpC0OmldtPPIy5sEBtS4EVohorT?usp=sharing))
                -   `complete`
                -   `hypothetical`
                -   `nucleotide`
                -   `rbp`
            -   `16Sep2022_data_excluding_refseq.tsv`
            -   `16Sep2022_phages_downloaded_from_genbank.gb` ([Download](https://drive.google.com/file/d/14LG1iGa1CqPbAjofZT1EY8VKnE8Iy45Q/view?usp=sharing))
        -   `models` ([Download](https://drive.google.com/drive/folders/1U5ugmkhD4LHElYnLj3B8Xt2TcPx-TOjB?usp=sharing))
            -   `boeckaerts.joblib`
            -   `esm.joblib`
            -   ...
        -   `preprocessing`
        -   `rbp_prediction`
        -   `temp`
        -   `1. Sequence Preprocessing.ipynb`
        -   ...
        -   `ClassificationUtil.py`
        -   ...

</details>

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

### Dependencies

**Operating System**: Windows, Linux, or macOS

Create a virtual environment with all the necessary dependencies installed via Conda (we recommend using [Miniconda](https://docs.anaconda.com/free/miniconda/index.html)):

```
conda env create -f environment_experiments.yaml
```

Activate this environment by running:

```
conda activate PHIEmbed-experiments
```

_Thanks to Dr. Paul K. Yu for sharing his environment configuration._

<details>
    <summary>Click here to show/hide note on running the notebook for protein embedding generation.
    </summary>

<br>

The notebook [`4. Protein Embedding Generation.ipynb`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/experiments/4.%20Protein%20Embedding%20Generation.ipynb) has a dependency (`bio_embeddings`) that requires it to be run on Unix or a Unix-like operating system. If you are using Windows, consider using [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) (WSL) or a virtual machine. We did not include `bio_embeddings` in [`environment_experiments.yaml`](https://github.com/bioinfodlsu/phage-host-prediction/blob/main/environment_experiments.yaml) to maintain cross-platform compatibility; you have to install it following the instructions [here](https://docs.bioembeddings.com/v0.2.3/).

Moreover, generating protein embeddings should ideally be done on a machine with a GPU. The largest (and best-performing) protein language model that we used, ProtT5, consumes 5.9 GB of GPU memory. If your local machine does not have a GPU or if its GPU has insufficient memory, we recommend using a cloud GPU platform.

**UPDATE (12 Jun 2023)**: In May 2023, Google Colab upgraded its Python runtime, resulting in compatibility issues with `bio_embeddings`. An alternative cloud GPU platform is Paperspace, which provides a [PyTorch 1.12 runtime](https://docs.digitalocean.com/products/paperspace/notebooks/concepts/notebook-templates/) that is compatible with `bio_embeddings`.

</details>

<details>
  <summary>Click here to show/hide the complete list of Python libraries and modules used in this project (excluding those that are part of the Python Standard Library).</summary>

<br>

| Library/Module                                                   | Description                                                                                                                                                                                | License                                           |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------- |
| [`pyyaml`](https://pyyaml.org/)                                  | Supports standard YAML tags and provides Python-specific tags that allow to represent an arbitrary Python object                                                                           | MIT License                                       |
| [`jsonnet`](https://jsonnet.org/)                                | Domain-specific language for JSON                                                                                                                                                          | Apache License 2.0                                |
| [`protobuf`](https://protobuf.dev/)                              | Google's language-neutral, platform-neutral, extensible mechanism for serializing structured data                                                                                          | BSD 3-Clause "New" or "Revised" License           |
| [`regex`](https://github.com/mrabarnett/mrab-regex)              | Provides additional functionality over the standard [`re`](https://docs.python.org/3/library/re.html) module while maintaining backwards-compatibility                                     | Apache License 2.0                                |
| [`nltk`](https://www.nltk.org/)                                  | Provides interfaces to corpora and lexical resources, along with a suite of text processing libraries for classification, tokenization, stemming, tagging, parsing, and semantic reasoning | Apache License 2.0                                |
| [`biopython`](https://biopython.org/)                            | Provides tools for computational molecular biology                                                                                                                                         | Biopython License Agreement, BSD 3-Clause License |
| [`ete3`](http://etetoolkit.org/)                                 | Provides functions for automated manipulation, analysis, and visualization of phylogenetic trees                                                                                           | GNU General Public License v3.0                   |
| [`pandas`](https://pandas.pydata.org/)                           | Provides functions for data analysis and manipulation                                                                                                                                      | BSD 3-Clause "New" or "Revised" License           |
| [`numpy`](https://numpy.org/)                                    | Provides a multidimensional array object, various derived objects, and an assortment of routines for fast operations on arrays                                                             | BSD 3-Clause "New" or "Revised" License           |
| [`scipy`](https://scipy.org/)                                    | Provides efficient numerical routines, such as those for numerical integration, interpolation, optimization, linear algebra, and statistics                                                | BSD 3-Clause "New" or "Revised" License           |
| [`scikit-learn`](https://scikit-learn.org/)                      | Provides efficient tools for predictive data analysis                                                                                                                                      | BSD 3-Clause "New" or "Revised" License           |
| [`xgboost`](https://xgboost.readthedocs.io/en/stable/)           | Implements machine learning algorithms under the gradient boosting framework                                                                                                               | Apache License 2.0                                |
| [`imbalanced-learn`](https://imbalanced-learn.org/stable/)       | Provides tools when dealing with classification with imbalanced classes                                                                                                                    | MIT License                                       |
| [`joblib`](https://joblib.readthedocs.io/en/latest/)             | Provides tools for lightweight pipelining in Python                                                                                                                                        | BSD 3-Clause "New" or "Revised" License           |
| [`cudatoolkit`](https://developer.nvidia.com/cuda-toolkit)       | Parallel computing platform and programming model for general computing on GPUs                                                                                                            | NVIDIA Software License                           |
| [`bio_embeddings`](https://docs.bioembeddings.com/v0.2.3/)       | Provides an interface for the use of language model-based biological sequence representations for transfer-learning                                                                        | MIT License                                       |
| [`torch`](https://pytorch.org/)                                  | Optimized tensor library for deep learning using GPUs and CPUs                                                                                                                             | BSD 3-Clause "New" or "Revised" License           |
| [`transformers`](https://huggingface.co/docs/transformers/index) | Provides pretrained models to perform tasks on different modalities such as text, vision, and audio                                                                                        | Apache License 2.0                                |
| [`sentencepiece`](https://github.com/google/sentencepiece)       | Unsupervised text tokenizer and detokenizer mainly for neural network-based text generation systems                                                                                        | Apache License 2.0                                |
| [`matplotlib`](https://matplotlib.org/)                          | Provides functions for creating static, animated, and interactive visualizations                                                                                                           | Matplotlib License (BSD-Compatible)               |
| [`umap-learn`](https://umap-learn.readthedocs.io/en/latest/)     | Implements uniform manifold approximation and projection, a dimensionality reduction technique                                                                                             | BSD 3-Clause "New" or "Revised" License           |

_The descriptions are taken from their respective websites._

</details>

↑ _Return to [Table of Contents](https://github.com/bioinfodlsu/phage-host-prediction?tab=readme-ov-file#table-of-contents)._

## 💻 Authors

-   **Mark Edward M. Gonzales** <br>
    gonzales.markedward@gmail.com

-   **Ms. Jennifer C. Ureta** <br>
    jennifer.ureta@gmail.com
    
-   **Dr. Anish M.S. Shrestha** <br>
    anish.shrestha@dlsu.edu.ph

This is a research project under the [Bioinformatics Laboratory](https://bioinfodlsu.com/), [Advanced Research Institute for Informatics, Computing and Networking](https://www.dlsu.edu.ph/research/research-centers/adric/), De La Salle University, Philippines.

This research was partly funded by the [Department of Science and Technology &ndash; Philippine Council for Health Research and Development](https://www.pchrd.dost.gov.ph/) (DOST-PCHRD) under the [e-Asia JRP 2021 Alternative therapeutics to tackle AMR pathogens (ATTACK-AMR) program](https://www.the-easia.org/jrp/projects/project_76.html). The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

[badge-jupyter]: https://img.shields.io/badge/Jupyter-F37626.svg?&style=flat&logo=Jupyter&logoColor=white
[badge-python]: https://img.shields.io/badge/python-3670A0?style=flat&logo=python&logoColor=white
[badge-pandas]: https://img.shields.io/badge/Pandas-2C2D72?style=flat&logo=pandas&logoColor=white
[badge-numpy]: https://img.shields.io/badge/Numpy-777BB4?style=flat&logo=numpy&logoColor=white
[badge-scipy]: https://img.shields.io/badge/SciPy-654FF0?style=flat&logo=SciPy&logoColor=white
[badge-github-actions]: https://img.shields.io/badge/GitHub_Actions-2088FF?style=flat&logo=github-actions&logoColor=white

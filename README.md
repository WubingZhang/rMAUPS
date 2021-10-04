## Introduction

Advances in proteomic profiling have enabled the study of protein regulation and degradation in diseases such as cancer. However, there are major computational challenges, such as how to perform quality control (QC) and normalization, how to identify differential protein abundance at multiple-levels of resolution (e.g. protein complexes), and how to integrate data with other omic technologies. Here, we developed a computational analysis pipeline Model-based Analysis of the Ubiquitin-Proteasome System using R (**rMAUPS**), which performs computational analysis of proteomics data efficiently and effectively. rMAUPS comprises four significant modules, including quality control (QC), imputation, differential analysis, and integrative analysis.

## Installation
Installing the package in a fresh environment may take a long time. It may fail because of some issues. You can check the possible issues and solutions from https://github.com/WubingZhang/rMAUPS/issues/3, or post a new issue there.

### Prerequisites
To install rMAUPS, you have to first install conda following the document (https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html#install-macos-silent).  The fast way is shown as follows:  

```
# Installing conda on macOS
$ wget wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
$ bash Miniconda3-latest-MacOSX-x86_64.sh

# Installing conda on Linux
$ wget wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh

# Installing conda on Windows
# 1. Download the installer (.exe file) from https://conda.io/miniconda.html
# 2. Install miniconda by double clicking the file
# 3. Test your installation. In your terminal window or Anaconda Prompt, run the command conda list.
```

You should also install r, r-recommended. If you have R installed before, you need to ensure libgfortran, libnetcdf and libxml2 are also installed in your conda environment. Besides, it seems important to have r-data.table and r-rcpparmadillo installed through conda before biocmanager installs dependencies (such as DESeq2).
```
$ conda install -c r r r-markdown r-recommended
$ conda install -c anaconda libnetcdf libxml2
$ conda install -c conda-forge r-devtools r-biocmanager r-ggpubr r-metap r-data.table r-rcpparmadillo
$ conda install -c bioconda bioconductor-gsva bioconductor-biocstyle
```

### Installation of rMAUPS in R
```
> install.packages(c("devtools", "BiocManager"), repos = "https://cloud.r-project.org")
# Install dependencies
> BiocManager::install(c("ggrepel", "DESeq2", "limma", "impute", "biomaRt", "msigdbr", "msmsTests"))
# Install rMAUPS from github
> devtools::install_github("WubingZhang/rMAUPS")
```

## Tutorial
Browse the detail documentation by running the command line below.
You can also visit http://cistrome.org/~wubing/rMAUPS.html and http://cistrome.org/~wubing/ubiquitin.html
```
> vignette("rMAUPS", package = "rMAUPS")
> vignette("ubiquitin", package = "rMAUPS")
```


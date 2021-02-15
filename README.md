# BowfinGOPhylogeny

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/eparey/binder_test/HEAD)

BowfinGOPhylogeny is a repository to build distance-based gene-order phylogenies. The notebook shows how to use the python3 scripts stored in `modules/` in conjonction with standard tree building tools to reproduce the phylogeny of 13 ray-finned fish presented in the [the bowfin genome paper](https://www.researchsquare.com/article/rs-92055/v1).

## Quick start with binder

**Click the binder badge above to directly open and interactively run the notebook.**

## Running locally

This is a short guide to setup the notebook environment with conda. If you do not wish to use conda you will need to setup the environment with the dependencies listed in `binder/environment.yaml`.

- Clone the repository

  ```
  git clone 
  cd 
  ```

- Install jupyter

  ```
  ```

- Create the conda environment

  We recommend using [Mamba](https://github.com/mamba-org/mamba) for a faster installation:

  ```
  conda install -c conda-forge mamba
  mamba env create -f binder/environment.yaml
  ```

  **Alternatively,** you can use conda directly :

  ```
  conda env create -f binder/environment.yaml
  ```

- Add the conda environment to jupyter kernels

  Install the extension
  ```
  ```

  Add

  ```
  ```

- Launch the notebook

  ```
  ```

## Authors

* [**Elise Parey**](mailto:elise.parey@bio.ens.psl.eu)
* **Ingo Braasch**
* **Hugues Roest Crollius**
* **Camille Berthelot**

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS.

    - [LICENSE-GPLv3.txt](LICENSE-GPLv3.txt)
    - [LICENSE-CeCILL.txt](LICENSE-CeCILLv2.txt)

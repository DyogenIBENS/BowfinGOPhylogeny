# BowfinGOPhylogeny

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/eparey/binder_test/HEAD)

**BowfinGOPhylogeny** is a repository to build distance-based gene-order phylogenies, based on the approach we used in [the bowfin genome paper](https://www.researchsquare.com/article/rs-92055/v1). The notebook shows how to use the python3 scripts stored in `modules/` in conjonction with standard tree building tools to reproduce the 13 ray-finned fish phylogeny we presented.

## Quick start with binder

Click the **launch binder** badge above to directly open and interactively run the notebook.

## Running locally

This is a short guide to setting up the notebook environment with conda. If you do not wish to use conda you will need to install the dependencies listed in `binder/environment.yaml` to run the notebook.

- Clone the repository and go to the root folder

```
git clone https://github.com/DyogenIBENS/BowfinGOPhylogeny.git
cd BowfinGOPhylogeny
```

- Create the conda environment using mamba

```
conda install -c conda-forge mamba
mamba env create -f binder/environment.yaml
```

- Install the conda extension for jupyter in your base environment - or in any environment where you want to run jupyter (if you do not have jupyter you should install it first with `mamba install -c anaconda jupyter`).

```
mamba install nb_conda_kernels
```

- Launch the notebook and then interactively select the bowfin_go kernel in the list

```
jupyter notebook
```

## Authors

* [**Elise Parey**](mailto:elise.parey@bio.ens.psl.eu)
* **Ingo Braasch**
* **Hugues Roest Crollius**
* **Camille Berthelot**

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS:

- [LICENSE-GPLv3.txt](LICENSE-GPLv3.txt)
- [LICENSE-CeCILL.txt](LICENSE-CeCILLv2.txt)

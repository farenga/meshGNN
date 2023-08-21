# DL-based Methods for Polygonal Grids Agglomeration

This repository contains the material related to the implementation of the GNN-based method introduced in the paper: [Agglomeration of Polygonal Grids using Graph Neural Networks with applications to Multigrid solvers](https://arxiv.org/abs/2210.17457).

P. F. Antonietti, N. Farenga, E. Manuzzi, G. Martinelli, L. Saverio

## Installation

To install python requirements run the command
    
    pip install -r requirements.txt

However, a dedicated installation of the torch-geometric package is recommended, via

    pip install torch-scatter -f https://data.pyg.org/whl/torch-{TORCH}+{CUDA}.html
    pip install torch-sparse -f https://data.pyg.org/whl/torch-{TORCH}+{CUDA}.html
    pip install torch-cluster -f https://data.pyg.org/whl/torch-{TORCH}+{CUDA}.html
    pip install torch-spline-conv -f https://data.pyg.org/whl/torch-{TORCH}+{CUDA}.html
    pip install torch-geometric

By replacing ``{TORCH}`` & ``{CUDA}`` with the installed ``torch`` and ``cuda`` versions, 
retrievable by running the python command ``torch.__version__``.

## Project structure

The folder structure is the following:
- [`dataset_generation`](dataset_generation) contains the MATLAB scripts used for generating the meshes, extracting the graphs together with their features, and the concatenation scripts to bundle the data.
- [`training`](training) contains the python scripts used for the definition, tuning and training of the models, the trained models, and the runtime benchmarking test.
- [`code_agglom_MG`](code_agglom_MG) contains the MATLAB code to perform agglomeration and numerical experiments (multigrid), with the different python wrappers to load the models within the `code_agglom_MG/mesh/agglomerate` folder. 

    > **Note**
    > To load a specific model, the path to the trained model at line 7 of the file `code_agglom_MG/mesh/agglomerate/aggl_GNN_fun.m`, and the model's wrapper filename at line 21, have to be both updated.

## Citation
    
    @misc{antonietti2022agglomeration,
        title={Agglomeration of Polygonal Grids using Graph Neural Networks with applications to Multigrid solvers},
        author={P. F. Antonietti and N. Farenga and E. Manuzzi and G. Martinelli and L. Saverio},
        year={2022},
        eprint={2210.17457},
        archivePrefix={arXiv},
        primaryClass={math.NA}
        }
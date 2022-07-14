# Interlimb coordination in Parkinson’s Disease is affected by a visuospatial dual task

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6835767.svg)](https://doi.org/10.5281/zenodo.6835767)

This results of this research are fully reproducible using the source code, Jupyter
notebooks, and raw data used to produce the results for the above titled project.

Raw data included in this repository:

- raw motion capture trials, as `.c3d` files
- system output from CAREN control software, D-Flow (includes treadmill speed, etc)
- OpenSim models

## Instructions

To run this analysis on your computer, both Julia and Jupyter (notebook or lab) must be available. A
version of Julia appropriate for your OS can be downloaded from the [Julia language
website](https://julialang.org/downloads/), and Jupyter can be installed from within Julia
(in the REPL) with

```julia
] add IJulia
```

Alternate instructions for installing Jupyter can be found on the [IJulia
github](https://github.com/JuliaLang/IJulia.jl) or the [Jupyter
homepage](https://jupyter.org/install) (advanced).

From within the main repository directory, start Julia and then start Jupyter in the Julia
REPL

```julia
using IJulia
notebook(;dir=pwd())
```

or if using a system Jupyter installation, start Jupyter from your favorite available shell
(e.g. Powershell on Windows, bash on any *nix variant, etc.).

The primary analysis is found in the `Analysis` notebook. Two other notebooks related to
this research are also included.


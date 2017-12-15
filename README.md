# BST 235 Final Project: Link Violation in GLMs

Andy Shi, Linglin Huang, Helian Feng

Final project for BIOSTAT 235, Fall 2017, Harvard T.H. Chan School of
Public Health.

## How to Install

This is an R package developed for the simulation part of our BST 235
Project. You can install the package using

    if (!require(devtools)) {
        install.packages("devtools")
    }
    devtools::install_github("shiandy/bst235Project")

## Guide

The main simulation functions are found in `R/simulation.R`, with some
helper functions in `R/fit_true_model.R`. The code to recreate the
figures is in the `vignettes` folder. To recreate the figures, first you
need to run `run-simulation.R`, which will run the simulation and store
the results in `data/simdata.csv`. The simulation is configured to run
on 4 cores/CPUs in parallel, and should take about 2.5 hours to
complete. Then, running the vignette in `bst235_finalproject.Rmd` will
generate all the necessary figures/plots.  The vignette also runs a
small experiment, so running the code in the vignette may also take some
time.

# Notebooks

This directory contains some Jupyter notebooks for one-off analyses or computational tools. Check below for documentation on each.

## TrainEBMVariantAnnotations

This notebook was created by Eyad Alqaysi as an internship project in Summer 2024. The goal of the notebook was to train EBM machine learning models to predict which reference-based annotations of variants were the strongest indicators of making false-positive or false-negative mistakes in variant calling, and comparing across some sequencing technologies. Users can use the `AnnotateVCF` WDL under Utilities to create comprehensive labels on variants that can be used to train these models. The notebook included here was specific to some data used during the summer, but can serve as inspiration for similar analyses.
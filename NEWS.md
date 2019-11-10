LUCIDus Updates
================

# LUCIDus 1.0.0

This is a feature update expanding the integrative clustering introduced in 0.9.0 to handle missing values in biomarker data as well as including an approach to estimate the standard errors of parameter estimates through the bootstrap method.

## New features

* Updated the `est_lucid()` function to handle missing values in biomarker data;
* Added the `boot_lucid()` function to estimate SE through bootstrapping.

## Other changes

* Added a testing incomplete data with missing biomarkers in the package and examples in the `est_lucid()` and `boot_lucid()` documentation;
* Updated the citation after getting published by the *[Bioinformatics](https://doi.org/10.1093/bioinformatics/btz667)*;
* Minor bug fixes.

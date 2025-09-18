# RDynLib_ftms

Semi-automated generation of reference spectral libraries from FTMS data

This folder contains the following files:

- [centroiding.qmd](centroiding.qmd) : in this file we convert the data from in
  profile mode to centroided mode.

- [ftms_preprocessing.qmd](ftms_preprocessing.qmd) : where we load and analyse
  the data with the *xcms* package.

- [ftms_filtering.qmd](ftms_filtering.qmd): in this file we filter the ftms data
  to keep one tree per feature.

- [ftms_assembled.qmd](https://github.com/rebdau/RDynLib_ftms/blob/ahlam/ftms_assembled.qmd) : In this quarto document
 we create an assembled MS/MS data from ftms object that contains one tree
 per feature created in the ftms_filtering.qmd file.
  
**The execution order:**

1. *centroiding.qmd*
2. *ftms_preprocessing.qmd*
3. *ftms_filtering.qmd*


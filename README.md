# RDynLib_ftms

Semi-automated generation of reference spectral libraries from FTMS data

This folder contains the following files:

- [centroiding.qmd](centroiding.qmd) : in this file we convert the data from in
  profile mode to centroided mode.

- [ftms_preprocessing.qmd](ftms_preprocessing.qmd) : where we load and analyse
  the data with the *xcms* package.

- [ftms_filtering.qmd](ftms_filtering.qmd): in this file we filter the ftms data
  to keep one tree per feature.

- [Add_new_mzml.qmd](https://github.com/rebdau/RDynLib_ftms/blob/ahlam/Add_new_mzml.qmd) : this workflow demonstrate how to process and add 
  a single mzml file to an existing sql database.

- [ftmsneg_flax_sql.qmd](https://github.com/rebdau/RDynLib_ftms/blob/main/ftmsneg_flax_sql.qmd) : 
 In this Quarto document, we create an SQL database of flax FTMS negative mode 
 data using the objects produced in ftms_filtering.qmd and ftms_assembled.qmd. 
 The database contains two types of spectra: assembled and separated.
 
**The execution order:**

1. *centroiding.qmd*
2. *ftms_preprocessing.qmd*
3. *ftms_filtering.qmd*
4. *ftmsneg_flax_sql.qmd*
5. *Add_new_mzml.qmd*
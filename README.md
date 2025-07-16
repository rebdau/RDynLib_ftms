# RDynLib_ftms

Semi-automated generation of reference spectral libraries from ftms data

This folder contains the following files:

-   " centroiding.qmd" : in this file we convert the data from in profile mode to centroided mode.

<!-- -->

-   " ftms_preprocessing.qmd" : where we upload and analyse the data with the xcms treatment.

-   "ftms_filtering.qmd" : in this file we filtered the ftms data and we kept just one tree per feature and we load the resulting object as "ftms_one_tree".

    **The execution order:**

    1.  " centroiding.qmd"
    2.  " ftms_preprocessing.qmd"
    3.  "ftms_filtering.qmd"

# Morgani _et al._ (2018) _Developmental Biology_

This is the repository associated with the article by Morgani _et al._, _A Sprouty4 reporter to monitor FGF/ERK signaling activity in ESC and mice_ published in **Developmental Biology** 441:1 pp104:126 (2018) [doi.org/10.1016/j.ydbio.2018.06.017](https://doi.org/10.1016/j.ydbio.2018.06.017)

This repository contains the corrected data tables for all preimplantation embryos analyzed in the study, the R scripts used to clean, transform and analyze the data, as well as the code to generate all plots in the article. All raw data tables, corrected data tables (also), original microscopy images and segmentation files can be found in Figshare [link to come]

This dataset contains data collected by both myself and Vidur Garg. File sets are labelled with the respective initials. Because of minor variation in data due to the experimenter, these data are processed separately (although in equivalent manner). Several scripts are thus duplicate, with _vg versions corresponding to Vidur Garg's data. Nestor Saiz contributed data in Figures 3, 5 and associated Figures S2 and S5; Vidur Garg contributed data in Figure 4 and associated Figure S4.

Scripts should work out of the box when using the .Rproj file and setting the working directory to the relevant folder. Please contact me with any issues or if any file is missing.

## Files included

**Folders**

* cor_files_ns: corrected data tables for all individual embryos analyzed by Nestor Saiz.
* cor_files_vg: corrected data tables for all individual embryos analyzed by Vidur Garg.
* movie_data: data tables with fluorescence measurements from timelapse movies analyzed in the study. 

**Data tables**

* spry4_exp_ref.csv: metadata table with experimental information.
* spry4_if.csv: metadata table with immunofluorescence details.
* spry4_mov_exp_ref.csv: metadata table with experimental information for live imaging experiments.
* spry-all-processed.csv: data table with all files combined, cleaned and transformed - output of ```read_data.R``` 

**Scripts**

* read_data.R: script to load corrected immunofluorescence data tables (from /cor_files), bind them into a single table, clean up, perform basic analyses, corrections and cell counts. 
* identify_spry.R: script to assign lineage identity to ICM cells using a basic k-means clustering approach, analogous to that we used in [Saiz _et al._ (2016)](https://www.nature.com/articles/ncomms13463) (see [saiz-et-al_2016](https://github.com/nestorsaiz/saiz-et-al_2016) repo).
* read_moviedata.R: script to load, do basic transformations and plot data from live imaging experiments (/movie_data). 
* do_counts: custom function to calculate total cell number, ICM cell number and litter average cell number for the given dataset.
* eb_cor.R: custom function to correct Z-associated fluorescence decay by fitting a linear regression to the log-transformed values and further correcting using an Emprirical Bayes approach, as we did in Saiz _et al._ (2016) (link in saiz-et-al_2016 repo README.md, see Methods there and script annotations).
* modas.R: custom function to calculate the mode for fluorescence Channel 2 (Venus or anti-GFP) for ICM cells.
* stage.R: custom function to assign embryos to developmental stage categories based on their total cell count.
* plotting-aes.R: contains objects used for aesthetics across figures.
* figure_3.R, figure_4.R and figure_5.R: scripts to generate the plots in the corresponding figures and associated supplementary figures.


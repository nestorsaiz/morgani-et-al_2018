# Morgani _et al._ (2018) _Developmental Biology_

This is the README file for the morgani-et-al_2018 repository associated with the article by Morgani _et al._, _A Sprouty4 reporter to monitor FGF/ERK signaling activity in ESC and mice_ **Developmental Biology** _(in press)_ [insert link here]

This repository contains the R scripts used to clean, transform and analyze the data, as well as the code to generate all plots in the article. All raw data tables, corrected data tables, original microscopy images and segmentation files can be found in Figshare [insert link]

## Files included

Data tables

* spry4_exp_ref.csv: metadata table with experimental information.
* spry4_if.csv: metadata table with immunofluorescence details.
* spry4_mov_exp_ref.csv: metadata table with experimental information for live imaging experiments.

Scripts

* read_data.R: script to load corrected immunofluorescence data tables (from /cor_files), bind them into a single table, clean up, perform basic analyses, corrections and cell counts. 
* identify_spry.R: script to assign lineage identity to ICM cells using a basic k-means clustering approach, analogous to that we used in Saiz _et al._ (2016) (see saiz-et-al_2016 repo).
* read_moviedata.R: script to load, do basic transformations and plot data from live imaging experiments (/movie_data). 
* do_counts: custom function to calculate total cell number, ICM cell number and litter average cell number for the given dataset.
* eb_cor.R: custom function to correct Z-associated fluorescence decay by fitting a linear regression to the log-transformed values and further correcting using an Emprirical Bayes approach, as we did in Saiz _et al._ (2016) (link in saiz-et-al_2016 repo README.md, see Methods there and script annotations).
* modas.R: custom function to calculate the mode for fluorescence Channel 2 (Venus or anti-GFP) for ICM cells.
* stage.R: custom function to assign embryos to developmental stage categories based on their total cell count.
* plotting-aes.R: contains objects used for aesthetics across figures.
* figure_3.R and figure_5.R: scripts to generate the plots in the corresponding figures and associated supplementary figures.

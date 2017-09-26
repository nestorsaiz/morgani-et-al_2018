## This function performs three basic calculations on the raw dataset:
## 1. Total cell count per embryo
## 2. Number of ICM cells per embryo
## 3. Average cell count per litter

library("plyr")
library("dplyr")

do.counts <- function(dataset) {
        ## Count the number of cells per embryo and add to main table
        counts <- dataset %>% group_by(Embryo_ID) %>% 
                summarize(Cellcount = n())
        dataset <- merge(dataset, counts)
        ## Calculate the number of ICM cells per embryo
        icmcounts <- dataset %>% filter(TE_ICM %in% c('ICM', 'in')) %>% 
                group_by(Embryo_ID) %>% 
                summarize(icm.count = n())
        dataset <- merge(dataset, icmcounts)
        ## Calculate the average cellcount per litter
        avg.litter <- dataset %>% group_by(Litter, Treatment) %>% 
                summarize(litter.mean = mean(Cellcount))
        ## Combine with main table and remove avg table
        dataset <- merge(dataset, avg.litter)
        return(dataset)
}
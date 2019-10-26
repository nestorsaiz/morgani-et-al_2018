## Load necessary packages
library('plyr')
library('dplyr')
library('reshape2')
library('ggplot2')
library('locfit')

################################################
## Read in data and calculate cell numbers #####
################################################

## Read in corrected .csv files
setwd("~/Documents/Data_Spry4/spry4_analysis/cor_files_vg")
files <- dir()
# Bind them all into a single file
spry.vg <- do.call(rbind.fill, lapply(files, read.csv))
rm(files)
# Reset working directory
setwd("~/Documents/Data_Spry4/spry4_analysis")
# Drop empty columns
spry.vg$X.1 <- as.character(spry.vg$X.1)
spry.vg$X.1[is.na(spry.vg$X.1)] <- " "
spry.vg$Cell_cycle <- ifelse(spry.vg$X.1 == 'div', 'M', 'I')
#spry$Cell_cycle <- as.factor(spry$Cell_cycle)
spry.vg$X.1 <- NULL

## Read in experimental reference file
spy.ref.vg <- read.csv('spry4_exp_ref_vg.csv')
## and combine exp.ref with main table
spry.vg <-  merge(spry.vg, spy.ref.vg)

names(spry.vg)[names(spry.vg) == "Embryo.Id"] <- "Embryo_ID"

## Calculate cell count per embryo in spry
spy.cellcount <- spry.vg %>% group_by(Embryo_ID) %>%
        summarize(Cellcount = n())
## and add count to main table
spry.vg <- merge(spry.vg, spy.cellcount)
rm(spy.cellcount)

## Load staging function and stage spry embryos
source('stage.R')
spry.vg <- stage(spry.vg)

## Calculate the average number of cells per litter
avg.litter <- spry.vg %>% group_by(Litter) %>% 
        summarize(litter.mean = mean(Cellcount))
## Combine with main table and remove avg table
spry.vg <- merge(spry.vg, avg.litter)
rm(avg.litter)

##Change numerical assignment to "TE" & "ICM"
spry.vg$TE_ICM <- ifelse(spry.vg$TE.ICM %in% c('1', '3', '4'), 'ICM', 'TE')

names(spry.vg)[names(spry.vg) == "TE.ICM"] <- "Id_man"
spry.vg$Id_man <- ifelse(spry.vg$Id_man == '1', 'EPI',
                      ifelse(spry.vg$Id_man == '2', 'TE',
                             ifelse(spry.vg$Id_man == '3', 'PrE', 'DP')))

spry.vg$Inlier.Outlier <- NULL
spry.vg <- rename(spry.vg, Cell_ID = Cell.ID)

################################################
## Correct Z-associated fluorescence decay #####
## and assign lineage identity to ICM cells ####
################################################

## Correct for Z-associated fluorescence decay for all channels
source('eb_cor.R')
channels <- c('CH1.Avg', 'CH2.Avg', 'CH3.Avg', 'CH4.Avg', 'CH5.Avg')
ebLogCor <- matrix(0, nrow = length(spry.vg$Embryo_ID), 
                   ncol = length(channels),
                   dimnames = list(c(), channels))
for (c in channels) {
        ebLogCor[, c] <- ebcor(spry.vg, c)
}
ebLogCor <- data.frame(ebLogCor)
ebLogCor <- rename(ebLogCor, CH1.ebLogCor = CH1.Avg, 
                   CH2.ebLogCor = CH2.Avg, 
                   CH3.ebLogCor = CH3.Avg,
                   CH4.ebLogCor = CH4.Avg,
                   CH5.ebLogCor = CH5.Avg)

## Combine spry with the EB corrected values
spry.vg <- cbind(spry.vg, ebLogCor)
rm(ebLogCor)

## Assign lineage identity to ICM cells depending on their location 
## in the log[GATA6] vs log[NANOG] space using k-means clustering
## see annotation on 'identify_spry.R' for details
#source('identify_spry_vg.R')

################################################
## Read in information on immunofluorescence ###
################################################

## Load immunofluorescence (IF) reference sheet
spy.if.vg <- read.csv('spry4_if_vg.csv')
## Extract info for Channels 2, 3 and 5 
## and their associated marker (primary ab) and secondary antibody (ab.2)
spy.aa <- subset(spy.if.vg, Channel %in% c('CH3', 'CH4', 'CH5')) %>% 
        group_by(Experiment, Litter, Channel, Marker) %>%
        summarize()
spy.bb <- subset(spy.if.vg, Channel %in% c('CH3', 'CH4', 'CH5')) %>% 
        group_by(Experiment, Litter, Channel, ab.2) %>%
        summarize()
## Cast in wide format 
spy.aa <- dcast(spy.aa, Experiment ~ Channel, value.var = 'Marker')
spy.aa <- rename(spy.aa, CH3.marker = CH3, CH4.marker = CH4, CH5.marker = CH5)
spy.bb <- dcast(spy.bb, Experiment ~ Channel, value.var = 'ab.2')
spy.bb <- rename(spy.bb, CH3.ab2 = CH3, CH4.ab2 = CH4, CH5.ab2 = CH5)
## and combine marker and ab.2 data
spy.if.vg <- merge(spy.aa, spy.bb)
rm(spy.aa, spy.bb)

# Merge with main table to divide embryos into
# those stained with anti-GFP and those with endogenous Venus
spry.vg <- merge(spry.vg, spy.if.vg)

################################################
## Tidy up and order factors for plotting ######
################################################

## Order factors as desired
spry$Genotype <- factor(spry$Genotype, levels = c('wt', 'het', 'ko'))
spry$Treatment <- NULL
spry$Tt_length <- NULL
spry$Tt_stage <- NULL
spry$Inlier.Outlier <- NULL
spry$venus.gfp <- NULL
#spry$venus.gfp <- factor(spry$venus.gfp, levels = c('Venus', 'GFP.ck', 'no.ab'))
spry$ab.2 <- NULL
#spry$ab.2 <- factor(spry$ab.2, levels = c('no.ab', 'Hoechst', 'af488.ck', 
                                          #'af568.gt', 'af647.rb', 'af647.ck'))
spry$Identity.km <- factor(spry$Identity.km, levels = c('DN', 'EPI', 'DP', 
                                                        'PRE', 'TE', 'morula'))
spry$TE_ICM <- factor(spry$TE_ICM, levels = c('TE', 'ICM', 'in', 'out'))

################################################

## Calculate number of embryos for each group 
## (stage x treatment x genotype x IF)
n.embryos <- spry %>% group_by(Embryo_ID, Stage, Genotype) %>% 
        summarize() %>% 
        group_by(Stage, Genotype) %>% 
        summarize(N = n())
## Write out the N numbers to a .csv file for quick reference
#n.embryos2 <- dcast(subset(n.embryos, Treatment != 'neg.control'), 
                    #Genotype + Stage ~ interaction(Treatment, venus.gfp), 
                    #value.var = 'N')
write.csv(n.embryos, file = 'N_numbers.csv', row.names = FALSE)

## If there are no errors in the code, 'All good :)' should print to the console
print('All good :)')

# ## Calculate the ratio of PrE to EPI per embryo
# ## Re-cast table to wide format and filter for ICM cells only
# spry.ratiocounts <- dcast(spry.lincounts, Experiment + Litter + 
#                                  Embryo_ID + TE_ICM + Cellcount + 
#                                  Stage + Exp_date + Img_date +
#                                  Genotype + Treatment + 
#                                  venus.gfp + ab.2 ~ Identity.km, 
#                          value.var = 'count')
# spry.ratiocounts <- spry.ratiocounts %>% filter(TE_ICM == 'ICM')
# ## Add EPI and DN numbers
# aa <- subset(spry.ratiocounts, DN >= 1)
# aa$DN.EPI <- aa$DN + aa$EPI
# bb <- subset(spry.ratiocounts, is.na(DN) == TRUE)
# bb$DN.EPI <- bb$EPI
# spry.ratiocounts <- rbind(aa, bb)
# rm(aa, bb)
# ## Select rows where both PRE and EPI are != 0 
# ## (otherwise cannot calculate ratio)
# spry.ratiocounts <- subset(spry.ratiocounts, PRE >= 1 & DN.EPI >= 1)
# ## Calculate the ratio as PRE/EPI
# spry.ratiocounts$ratio <- spry.ratiocounts$PRE / spry.ratiocounts$DN.EPI
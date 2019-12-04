## Load necessary packages
library('plyr')
library('dplyr')
library('reshape2')
library('ggplot2')
library('locfit')

options(stringsAsFactors = F)

################################################
## Read in data and calculate cell numbers #####
################################################

## Read in Nestor Saiz's (NS) corrected .csv files
setwd("~/Documents/Data_Spry4/spry4_analysis/cor_files_ns")
files <- dir()
# Bind them all into a single file
spry.ns <- do.call(rbind.fill, lapply(files, read.csv))
rm(files)

## Read in Vidur Garg's (VG) corrected .csv files
setwd("~/Documents/Data_Spry4/spry4_analysis/cor_files_vg")
files <- dir()
# Bind them all into a single file
spry.vg <- do.call(rbind.fill, lapply(files, read.csv))
rm(files)

## Reset working directory
setwd("~/Documents/Data_Spry4/spry4_analysis")

## Extract cell cycle information in Vidur's files
spry.vg$X.1 <- as.character(spry.vg$X.1)
spry.vg$X.1[is.na(spry.vg$X.1)] <- " "
spry.vg$Cell_cycle <- ifelse(spry.vg$X.1 == 'div', 'M', 'I')

## Drop empty columns
spry.vg$X.1 <- NULL
spry.ns$X.1 <- NULL
spry.ns$X.2 <- NULL
spry.ns$X.3 <- NULL
spry.vg$Inlier.Outlier <- NULL

## Transform/rename variables from both datasets to match
# Create TE_ICM column based on manual identity assignment
spry.vg$TE_ICM <- ifelse(spry.vg$TE.ICM %in% c('1', '3', '4'), 'ICM', 'TE')
# Rename variables in Vidur's dataset
spry.vg <- rename(spry.vg, Embryo_ID = Embryo.Id, 
                  Id_man = TE.ICM, 
                  Cell_ID = Cell.ID, 
                  hoechst.Avg = CH1.Avg, 
                  bf.Avg = CH2.Avg, 
                  green.Avg = CH3.Avg, 
                  farred.Avg = CH4.Avg, 
                  red.Avg = CH5.Avg)
# Translate numerical manual identity assignment
spry.vg$Id_man <- ifelse(spry.vg$Id_man == '1', 'EPI',
                         ifelse(spry.vg$Id_man == '2', 'TE',
                                ifelse(spry.vg$Id_man == '3', 'PrE', 'DP')))
# Rename variables in Nestor's dataset
spry.ns <- rename(spry.ns, hoechst.Avg = CH1.Avg, 
                  green.Avg = CH2.Avg, 
                  farred.Avg = CH3.Avg, 
                  bf.Avg = CH4.Avg, 
                  red.Avg = CH5.Avg)

## Bind both datasets into one main table 
spry <- rbind.fill(spry.ns, spry.vg)

## Drop sum columns for all channels
spry$CH1.Sum <- NULL
spry$CH2.Sum <- NULL
spry$CH3.Sum <- NULL
spry$CH4.Sum <- NULL
spry$CH5.Sum <- NULL

## Read in experimental reference file
spy.ref <- read.csv('spry4_exp_ref.csv')
# and combine with main table
spry <-  merge(spry, spy.ref)

## Remove from dataset embryos with unknown genotypes
spry <- spry %>% filter(Genotype1 != 'unknown')
spry <- spry %>% filter(Genotype2 != 'unknown')

## Calculate total number of cells, 
## litter average and
## number of ICM cells per embryo
source('do_counts.R')
spry <- do.counts(spry)

## Calculate the average number of cells in het embryos
## for each litter (all litters have hets, but not all have wt)
## to compare size of different genotypes in each litter
het.meancount <- spry %>% filter(Experimenter == 'NS', 
                                Treatment == 'Littermate', 
                                Genotype1 == 'het') %>% 
        group_by(Litter) %>% 
        summarize(het.meancount = mean(Cellcount))

## Load staging function and stage spry embryos
source('stage.R')
spry <- stage(spry)

################################################
## Correct Z-associated fluorescence decay #####
## and assign lineage identity to ICM cells ####
################################################

## Correct for Z-associated fluorescence decay for all channels
source('eb_cor.R')
channels <- c('hoechst.Avg', 'green.Avg', 'farred.Avg', 'bf.Avg', 'red.Avg')
ebLogCor <- matrix(0, nrow = length(spry$Embryo_ID), 
                   ncol = length(channels),
                   dimnames = list(c(), channels))
for (c in channels) {
        ebLogCor[, c] <- ebcor(spry, c)
}
ebLogCor <- data.frame(ebLogCor)
ebLogCor <- rename(ebLogCor, hoechst.ebLogCor = hoechst.Avg, 
                   green.ebLogCor = green.Avg, 
                   farred.ebLogCor = farred.Avg,
                   bf.ebLogCor = bf.Avg,
                   red.ebLogCor = red.Avg)

## Combine spry with the EB corrected values
spry <- cbind(spry, ebLogCor)
rm(ebLogCor)

## Assign lineage identity to ICM cells depending on their location 
## in the log[GATA6] vs log[NANOG] space using k-means clustering
## see annotation on 'identify_spry.R' for details
source('identify_spry.R')

################################################
## Read in information on immunofluorescence ###
################################################

## Load immunofluorescence (IF) reference sheet
spy.if <- read.csv('spry4_if.csv')
## Extract info for Channels 2, 3 and 5 for Nestor's data 
## and their associated marker (primary ab) and secondary antibody (ab.2)
spy.aa <- subset(spy.if, Channel %in% c('CH2', 'CH3', 'CH5') &
                         Experimenter == 'NS') %>% 
        group_by(Experiment, Litter, Channel, Marker) %>%
        summarize()
spy.bb <- subset(spy.if, Channel %in% c('CH2', 'CH3', 'CH5') &
                         Experimenter == 'NS') %>% 
        group_by(Experiment, Litter, Channel, ab.2) %>%
        summarize()
## Cast in wide format 
spy.aa <- dcast(spy.aa, Experiment ~ Channel, value.var = 'Marker')
spy.aa <- rename(spy.aa, green.marker = CH2, 
                 farred.marker = CH3, 
                 red.marker = CH5)
spy.bb <- dcast(spy.bb, Experiment ~ Channel, value.var = 'ab.2')
spy.bb <- rename(spy.bb, green.ab2 = CH2, 
                 farred.ab2 = CH3, 
                 red.ab2 = CH5)
## and combine marker and ab.2 data
spy.if.ns <- merge(spy.aa, spy.bb)
rm(spy.aa, spy.bb)

## Extract info for Channels 3, 4 and 5 for Vidur's data 
## and their associated marker (primary ab) and secondary antibody (ab.2)
spy.aa <- subset(spy.if, Channel %in% c('CH3', 'CH4', 'CH5') &
                         Experimenter == 'VG') %>% 
        group_by(Experiment, Litter, Channel, Marker) %>%
        summarize()
spy.bb <- subset(spy.if, Channel %in% c('CH3', 'CH4', 'CH5') &
                         Experimenter == 'VG') %>% 
        group_by(Experiment, Litter, Channel, ab.2) %>%
        summarize()
## Cast in wide format 
spy.aa <- dcast(spy.aa, Experiment ~ Channel, value.var = 'Marker')
spy.aa <- rename(spy.aa, green.marker = CH3, 
                 farred.marker = CH4, 
                 red.marker = CH5)
spy.bb <- dcast(spy.bb, Experiment ~ Channel, value.var = 'ab.2')
spy.bb <- rename(spy.bb, green.ab2 = CH3, 
                 farred.ab2 = CH4, 
                 red.ab2 = CH5)
## and combine marker and ab.2 data
spy.if.vg <- merge(spy.aa, spy.bb)
rm(spy.aa, spy.bb)

## Merge modified spy.if for Nestor and Vidur into original
spy.if <- rbind(spy.if.ns, spy.if.vg)

# Merge with main table to divide embryos into
# those stained with anti-GFP and those with endogenous Venus
spry <- merge(spry, spy.if)

################################################
## Tidy up and order factors for plotting ######
################################################

## Order factors as desired
spry$Genotype1 <- factor(spry$Genotype1, levels = c('wt', 'het', 'homo', 
                                                    'unknown'))
spry$Genotype2 <- factor(spry$Genotype2, levels = c('wt', 'het', 'ko', 
                                                    'unknown'))
spry$Treatment <- factor(spry$Treatment, levels = c('Littermate', 'Control', 
                                                    'FGF4_1000', 'AZD_1', 
                                                    'PD03_1', 'neg.control'))
spry$green.marker <- factor(spry$green.marker, levels = c('Venus', 'GFP.ck', 
                                                          'no.ab'))
spry$green.ab2 <- factor(spry$green.ab2, levels = c('no.ab', 'af488.ck'))
spry$TE_ICM <- factor(spry$TE_ICM, levels = c('TE', 'ICM', 'out', 'in'))
spry$Tt_stage <- factor(spry$Tt_stage, levels = c('8cell', '64_90', '90_120'))

################################################

## Write main table out to file for other applications
write.csv(spry, file = 'spry-all-raw.csv', row.names = FALSE)

## Write out a table of 'littermates' (non-cultured embryos)
## stained for NANOG and GATA6 to use for data mining
lms <- subset(spry, Treatment == 'Littermate' & 
                      red.marker == 'GATA6.gt' & 
                      farred.marker == 'NANOG.rb' & 
                      Genotype2 == 'wt')
write.csv(lms, file = 'lms-spry-raw.csv', row.names = FALSE)

## Calculate number of embryos for each group 
## (stage x treatment x genotype x IF)
n.embryos <- spry %>% group_by(Embryo_ID, Stage, Treatment, Gene1, Genotype1, 
                               Gene2, Genotype2, green.marker) %>% 
        summarize() %>% 
        group_by(Stage, Treatment, Gene1, Genotype1, 
                 Gene2, Genotype2, green.marker) %>% 
        summarize(N = n())
## Write out the N numbers to a .csv file for quick reference
n.embryos2 <- dcast(subset(n.embryos, Treatment != 'neg.control'), 
                    Gene1 + Genotype1 + Gene2 + Genotype2 + Stage ~ 
                            interaction(Treatment, green.marker), 
                    value.var = 'N')
write.csv(n.embryos2, file = 'N_numbers.csv', row.names = FALSE)

## If there are no errors in the code, 'All good :)' should print to the console
print('All good :)')
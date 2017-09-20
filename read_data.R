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

## Read in corrected .csv files
setwd("~/Documents/Data_Spry4/spry4_analysis/cor_files")
files <- dir()
# Bind them all into a single file
spry <- do.call(rbind.fill, lapply(files, read.csv))
rm(files)
# Reset working directory
setwd("~/Documents/Data_Spry4/spry4_analysis")
# Drop empty columns
spry$X.1 <- NULL
spry$X.2 <- NULL
spry$X.3 <- NULL

## Read in experimental reference file
spy.ref <- read.csv('spry4_exp_ref.csv')
## Rename 'homo' genotype to 'ko' since it is a knockout
ko <- subset(spy.ref, Genotype1 == 'homo')
ko$Genotype1 <- 'ko'
rest <- subset(spy.ref, Genotype1 != 'homo')
rest$Genotype1 <- as.character(rest$Genotype1)
## Combine ko and rest to produce spy.ref again
spy.ref <- rbind(rest, ko)
## and combine exp.ref with main table
spry <-  merge(spry, spy.ref)

## Calculate total number of cells per embryo, litter average and
## number of ICM cells per embryo
source('do_counts.R')
spry <- do.counts(spry)

## Load staging function and stage spry embryos
source('stage.R')
spry <- stage(spry)

################################################
## Correct Z-associated fluorescence decay #####
## and assign lineage identity to ICM cells ####
################################################

## Correct for Z-associated fluorescence decay for all channels
source('eb_cor.R')
channels <- c('CH1.Avg', 'CH2.Avg', 'CH3.Avg', 'CH5.Avg')
ebLogCor <- matrix(0, nrow = length(spry$Embryo_ID), 
                   ncol = length(channels),
                   dimnames = list(c(), channels))
for (c in channels) {
        ebLogCor[, c] <- ebcor(spry, c)
}
ebLogCor <- data.frame(ebLogCor)
ebLogCor <- rename(ebLogCor, CH1.ebLogCor = CH1.Avg, 
                   CH2.ebLogCor = CH2.Avg, 
                   CH3.ebLogCor = CH3.Avg,
                   CH5.ebLogCor = CH5.Avg)

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
# Extract information for Channel 2 (Venus or anti-GFP)
# as rest of channels are the same across (except for neg.controls)
spy.gfp <- subset(spy.if, Channel == 'CH2')
# Collapse information to Experiment, litter and marker
spy.gfp <- spy.gfp %>% group_by(Experiment, 
                                Litter, 
                                Marker, 
                                ab.2) %>% 
        summarize()
spy.gfp <- rename(spy.gfp, CH2.marker = Marker, 
                  CH2.ab2 = ab.2)

# Merge with main table to divide embryos into
# those stained with anti-GFP and those with endogenous Venus
spry <- merge(spry, spy.gfp)

################################################
## Tidy up and order factors for plotting ######
################################################

## Order factors as desired
spry$Genotype1 <- factor(spry$Genotype1, levels = c('wt', 'het', 'ko', 'unknown'))
spry$Treatment <- factor(spry$Treatment, levels = c('Littermate', 'Control', 
                                                    'FGF4_1000', 'AZD_1', 
                                                    'PD03_1', 'neg.control'))
spry$CH2.marker <- factor(spry$CH2.marker, levels = c('Venus', 'GFP.ck', 'no.ab'))
spry$CH2.ab2 <- factor(spry$CH2.ab2, levels = c('no.ab', 'Hoechst', 'af488.ck', 
                                          'af568.gt', 'af647.rb', 'af647.ck'))
spry$TE_ICM <- factor(spry$TE_ICM, levels = c('TE', 'ICM', 'in', 'out'))

################################################

## Calculate number of embryos for each group 
## (stage x treatment x genotype x IF)
n.embryos <- spry %>% group_by(Embryo_ID, Stage, Treatment, 
                               Genotype1, CH2.marker) %>% 
        summarize() %>% 
        group_by(Stage, Treatment, Genotype1, CH2.marker) %>% 
        summarize(N = n())
## Write out the N numbers to a .csv file for quick reference
n.embryos2 <- dcast(subset(n.embryos, Treatment != 'neg.control'), 
                    Genotype1 + Stage ~ interaction(Treatment, CH2.marker), 
                    value.var = 'N')
write.csv(n.embryos2, file = 'N_numbers.csv', row.names = FALSE)

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
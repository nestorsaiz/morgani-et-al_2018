## Load necessary packages
library('plyr')
library('dplyr')
library('reshape2')
library('ggplot2')
library('locfit')

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

## Read in experimental reference file
spy.ref <- read.csv('spry4_exp_ref.csv')
# combine exp.ref with main table
spry <-  merge(spry, spy.ref)

## Calculate cell count per embryo in spry
spy.cellcount <- spry %>% group_by(Embryo_ID) %>%
        summarize(Cellcount = n())
spry <- merge(spry, spy.cellcount)
rm(spy.cellcount)

## Load staging function and stage spry
source('stage.R')
spry <- stage(spry)

## Calculate the average number of cells per litter
avg.litter <- spry %>% group_by(Litter) %>% 
        summarize(litter.mean = mean(Cellcount))
## Combine with main table and remove avg table
spry <- merge(spry, avg.litter)
rm(avg.litter)

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

## Assign identities using the thresholding method
source('identify.R')
spry <- id.linear(spry)

## Load immunofluorescence (IF) reference sheet
spy.if <- read.csv('spry4_if.csv')
# Extract information for Channel 2 (anti-GFP or Venus)
# as rest of channels are the same across (except for neg.controls)
spy.gfp <- subset(spy.if, Channel == 'CH2')
# Collapse information to Experiment, litter and marker
spy.gfp <- spy.gfp %>% group_by(Experiment, 
                                Litter, 
                                Marker, 
                                ab.2) %>% 
        summarize()
spy.gfp <- rename(spy.gfp, CH2.marker = Marker)

# Merge with main table to divide embryos into
# those stained with anti-GFP and those with endogenous Venus
spry <- merge(spry, spy.gfp)

## Order factors as desired
spry$Genotype <- factor(spry$Genotype, levels = c('wt', 'het', 'unknown', ''))
spry$Treatment <- factor(spry$Treatment, levels = c('Littermate', 'Control', 
                                                    'FGF4_1000', 'PD03_1', 
                                                    'AZD_1', 'neg.control'))
spry$CH2.marker <- factor(spry$CH2.marker, levels = c('Venus', 'GFP.ck'))
spry$ab.2 <- factor(spry$ab.2, levels = c('NA', 'Hoechst', 'af488.ck', 
                                          'af568.gt', 'af647.rb'))
spry$Treatment <- factor(spry$Treatment, levels = c('Littermate', 'Control', 
                                                    'FGF4_1000', 'AZD_1', 
                                                    'PD03_1', 'neg.control'))

## Calculate number of embryos for each group (stage x treatment x genotype x IF)
n.embryos <- spry %>% group_by(Embryo_ID, Stage, Treatment, 
                               Genotype, CH2.marker) %>% 
        summarize() %>% 
        group_by(Stage, Treatment, Genotype, CH2.marker) %>% 
        summarize(N = n())
## Write out the N numbers to a .csv file for quick reference
n.embryos2 <- dcast(subset(n.embryos, Treatment != 'neg.control'), 
                   Genotype + Stage ~ interaction(Treatment, CH2.marker), 
                   value.var = 'N')
write.csv(n.embryos2, file = 'N_numbers.csv', row.names = FALSE)

## Calculate the average level for each fluorescence channel
## for each embryo and lineage
meantensity <- spry %>%
        group_by(Experiment, Litter, Embryo_ID, Identity.lin, TE_ICM, 
                 Cellcount, Stage, Exp_date, Img_date, Genotype, 
                 Treatment, Tt_length, Tt_stage, CH2.marker, ab.2) %>%
        summarize(Count = n(), 
                  CH2.mean = mean(CH2.ebLogCor), 
                  CH3.mean = mean(CH3.ebLogCor), 
                  CH5.mean = mean(CH5.ebLogCor))

## Calculate the number of cells per lineage for each embryo
spry.lincounts <- spry %>% filter(Treatment == 'Littermate') %>% 
        group_by(Experiment, Litter, 
                 Embryo_ID, TE_ICM, 
                 Cellcount, Stage, 
                 Exp_date, Img_date, 
                 Genotype, Treatment, 
                 Identity.lin, CH2.marker, ab.2) %>%
        summarize(count = n())

## Calculate the ratio of PrE to EPI per embryo
## Re-cast table to wide format and filter for ICM cells only
spry.ratiocounts <- dcast(spry.lincounts, Experiment + Litter + 
                                 Embryo_ID + TE_ICM + Cellcount + 
                                 Stage + Exp_date + Img_date +
                                 Genotype + Treatment + 
                                 CH2.marker + ab.2 ~ Identity.lin, 
                         value.var = 'count')
spry.ratiocounts <- spry.ratiocounts %>% filter(TE_ICM == 'ICM')
## Add EPI and DN numbers
aa <- subset(spry.ratiocounts, DN >= 1)
aa$DN.EPI <- aa$DN + aa$EPI
bb <- subset(spry.ratiocounts, is.na(DN) == TRUE)
bb$DN.EPI <- bb$EPI
spry.ratiocounts <- rbind(aa, bb)
rm(aa, bb)
## Select rows where both PRE and EPI are != 0 
## (otherwise cannot calculate ratio)
spry.ratiocounts <- subset(spry.ratiocounts, PRE >= 1 & DN.EPI >= 1)
## Calculate the ratio as PRE/EPI
spry.ratiocounts$ratio <- spry.ratiocounts$PRE / spry.ratiocounts$DN.EPI

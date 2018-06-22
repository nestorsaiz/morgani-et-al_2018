## Load necessary packages
library('plyr')
library('dplyr')
library('reshape2')
library('ggplot2')
library('locfit')

options(stringsAsFactors = F)

################################################
## Read in data ################################
################################################

## Read in movie segmentation .csv files
setwd("~/Documents/Data_Spry4/spry4_analysis/movie_data")
files <- dir()
# Bind them all into a single file
spry.mov <- do.call(rbind.fill, lapply(files, read.csv))
rm(files)
## Reset working directory
setwd("~/Documents/Data_Spry4/spry4_analysis")

## Read in experimental reference file and merge with main table
mov.expref <- read.csv('spry4_mov_exp_ref.csv')
spry.mov <- merge(spry.mov, mov.expref)

## Remove data points where Area < 50 (often just noise or outliers)
spry.mov <- subset(spry.mov, Area > 50)

## Turn Frames (15 min timepoints) into hours
cuts <- seq(0, 72, by = 4)
hours <- seq(1, 18, by = 1)
spry.mov$hours <- cut(spry.mov$Frame, breaks = cuts, labels = hours)
spry.mov$hours <- as.numeric(spry.mov$hours)
rm(hours)

## Break timeline into 4h bins
cuts <- seq(0, 72, by = 16)
cuts <- c(cuts, Inf)
bins <- c('0-4h', '4-8h', '8-12h', '12-16h', '16+ h')
spry.mov$h4 <- cut(spry.mov$Frame, breaks = cuts, labels = bins)
rm(bins)

## Order treatment levels
spry.mov$Treatment <- factor(spry.mov$Treatment, 
                             levels = c('Control', 'AZD_1', 'PD03_1'))

## Calculate average intensity per hour and Z-slice, 
## across all controls in each litter
avg.control <- spry.mov %>% filter(Treatment == 'Control') %>% 
        group_by(Litter, hours, Slice) %>% 
        summarize(avg.control = mean(Mean))
## Combine with main data frame
spry.mov <- merge(spry.mov, avg.control)

avg.frame1 <- spry.mov %>% filter(Frame == 1) %>% 
        group_by(Litter, Treatment, Slice) %>% 
        summarize(avg.frame1 = mean(Mean))
spry.mov <- merge(spry.mov, avg.frame1)

## Drop inidivudal time frames and normalize the intensity value ('Mean') 
## of each Z-slice ('Slice'), hour and embryo against 
## the control average for that litter (avg.control)
spry.mov <- spry.mov %>% 
        group_by(Slice, Litter, hours, Embryo_ID, Area, Mean, Min, Max, 
                 Median, MinThr, MaxThr, Experiment, Treatment, h4, 
                 avg.frame1, Gene1, Genotype1) %>% 
        summarize()
spry.mov$meanorm <- spry.mov$Mean / spry.mov$avg.frame1

## Plot all intensity values (per slice, per frame, per embryo) in 4h bins

fig5f <- ggplot(data = spry.mov, 
                aes(x = h4, y = meanorm))
fig5f <- fig5f + geom_boxplot(color = 'black', 
                              outlier.shape = 1, outlier.size = 1)
fig5f <- fig5f + geom_jitter(aes(color = Slice), size = 0.75, 
                             position = position_jitter(width = 0.3), 
                             alpha = 0.5)
fig5f <- fig5f + looks + scale_color_distiller(direction = 1, palette = 'Blues')
fig5f <- fig5f + facet_wrap( ~ Treatment)
fig5f <- fig5f + labs(x = 'Hour bins', 
                      y = 'Normalized average intensity per z-slice, per hour')
print(fig5f)


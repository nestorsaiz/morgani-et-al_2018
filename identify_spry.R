## This script assigns lineage identity to ICM cells for the Spry4 dataset
## using a k-means clustering strategy 
## analogous to that we published in Saiz et al (2016) Nat Comms

## Function to do an exploratory plot of GATA6 vs NANOG values for ICM cells
spread.icm <- function(dataset) {
        ol <- ggplot(dataset, 
                     aes(x = red.ebLogCor, y = farred.ebLogCor))
        ol <- ol + geom_jitter(aes(color = Cellcount), size = 3, alpha = 0.6)
        ol <- ol + coord_fixed() + theme_bw()
        ol <- ol + theme(axis.text = element_text(size = 15), 
                         axis.title = element_text(size = 20), 
                         legend.text = element_text(size = 15), 
                         legend.title = element_text(size = 20))
        ol <- ol + scale_color_distiller(direction = 1, 
                                         palette = 'Blues')
        ol <- ol + labs(y = 'log[NANOG]', x = 'log[GATA6]')
        print(ol)
}

## Subset and combine ICM levels for Littermates between 32-64 and 90-150 cells
## The overlay of these populations generate discrete clusters 
## that can be used to find centers for DP, PrE and EPI cells
## Using the full range from 32-150 cells generates a continuous cloud
## where single clusters cannot be found

set.seed(20170401)

## Perform for Nestor's subset of data
sns <- subset(spry, Experimenter == 'NS')
bb <- sns %>% filter(Treatment != 'neg.control', 
                      Genotype2 == 'wt', 
                      ## Exclude litters not stained for GATA6 and NANOG
                      !Litter %in% c('AG', 'V'), 
                      ## Do not subset for wild type embryos
                      ## as data shows a similar distribution 
                      ## for wt and hets
                      ## more data will yield denser clusters
                      #Genotype ==  'wt', 
                      Stage %in% c('32_64', '90_120', '120_150'),  
                      TE_ICM %in% c('ICM', 'in')) %>% 
        select(Cellcount, red.ebLogCor, farred.ebLogCor)
## Scatter plot of selected data
spread.icm(bb)

## Calculate three centers for bb (DP, PrE and EPI)
pp <- kmeans(bb[, 2:3], 3)
## Extract the minimum center values for each channel 
## and append to the centers list as the DN center
dn <- c(pp$centers[, 1][which.min(pp$centers[, 1])], 
        pp$centers[, 2][which.min(pp$centers[, 2])])
centers.ns <- rbind(pp$centers, dn, deparse.level = 0)

## Create conditions to test for ICM (is.icm) and morulas (is.morula)
is.icm <- sns$Cellcount > 31 & sns$TE_ICM == 'ICM' & 
        !sns$Litter %in% c('V', 'AG')
is.te <- sns$Cellcount > 31 & sns$TE_ICM == 'TE'
is.morula <- sns$Cellcount <= 31

## Make matrix to hold sum of squares (rows = icm rows, 4 columns)
ssq <- matrix(0, length(sns$TE_ICM[is.icm]), 4)
## and populate with min sum of squares for each cell
## (min distance to each of the centers for each cell's GATA6 vs NANOG value)
for(i in 1:4) {
        ssq[,i] <- (sns$red.ebLogCor[is.icm] - centers.ns[i,1])^2 + 
                (sns$farred.ebLogCor[is.icm] - centers.ns[i,2])^2
}
## calculate what center each cell is closest to 
## (which sum of squares is smallest)
min.ssq <- apply(ssq, 1, which.min)

## Create new variable in sns to hold the lineage identity
## assigned using k-means clustering (Identity.km)
sns$Identity.km <- rep(NA, nrow(sns))
## Morulas and TE cells remain unchaged
sns$Identity.km[is.morula] <- 'morula'
sns$Identity.km[is.te] <- 'TE'
## Assign identity to ICM cells based on min.ssq values
sns$Identity.km[is.icm] <- c('DP', 'EPI', 'PRE', 'DN')[min.ssq]

## Perform for Vidur's subset of data
svg <- subset(spry, Experimenter == 'VG')
bb <- svg %>% filter(Treatment != 'neg.control', 
                     Genotype2 == 'wt', 
                     ## Do not subset for wild type embryos
                     ## as data shows a similar distribution 
                     ## for wt and hets
                     ## more data will yield denser clusters
                     #Genotype ==  'wt', 
                     Stage %in% c('32_64', '90_120', '120_150'),  
                     TE_ICM %in% c('ICM', 'in')) %>% 
        select(Cellcount, red.ebLogCor, farred.ebLogCor)
## Scatter plot of selected data
spread.icm(bb)

## Calculate three centers for bb (DP, PrE and EPI)
pp <- kmeans(bb[, 2:3], 3)
## Extract the minimum center values for each channel 
## and append to the centers list as the DN center
dn <- c(pp$centers[, 1][which.min(pp$centers[, 1])], 
        pp$centers[, 2][which.min(pp$centers[, 2])])
centers.vg <- rbind(pp$centers, dn, deparse.level = 0)

## Create conditions to test for ICM (is.icm) and morulas (is.morula)
is.icm <- svg$Cellcount > 31 & svg$TE_ICM == 'ICM'
is.te <- svg$Cellcount > 31 & svg$TE_ICM == 'TE'
is.morula <- svg$Cellcount <= 31

## Make matrix to hold sum of squares (rows = icm rows, 4 columns)
ssq <- matrix(0, length(svg$TE_ICM[is.icm]), 4)
## and populate with min sum of squares for each cell
## (min distance to each of the centers for each cell's GATA6 vs NANOG value)
for(i in 1:4) {
        ssq[,i] <- (svg$red.ebLogCor[is.icm] - centers.vg[i,1])^2 + 
                (svg$farred.ebLogCor[is.icm] - centers.vg[i,2])^2
}
## calculate what center each cell is closest to 
## (which sum of squares is smallest)
min.ssq <- apply(ssq, 1, which.min)

## Create new variable in svg to hold the lineage identity
## assigned using k-means clustering (Identity.km)
svg$Identity.km <- rep(NA, nrow(svg))
## Morulas and TE cells remain unchaged
svg$Identity.km[is.morula] <- 'morula'
svg$Identity.km[is.te] <- 'TE'
## Assign identity to ICM cells based on min.ssq values
svg$Identity.km[is.icm] <- c('DP', 'PRE', 'EPI', 'DN')[min.ssq]

## Combine both subsets into one again
spry <- rbind(sns, svg)

## Order factors as desired
spry$Identity.km <- factor(spry$Identity.km, levels = c('TE', 'PRE', 'DP', 
                                                        'EPI', 'DN', 'morula'))

## Load standard plotting aesthetics
source('plotting-aes.R')
## Plot data by stage to visualize the outcome
qplot(red.ebLogCor,  farred.ebLogCor,
      data = subset(spry, TE_ICM == 'ICM' & Treatment == "Littermate" & 
                            !Litter %in% c('AG', 'V') & Genotype1 != 'unknown' & 
                            Genotype2  == 'wt'), 
      color = Identity.km) + looks + scale_color_manual(values = idcols) + 
        facet_grid(Genotype1 + Experimenter ~ Stage) + coord_fixed()
## This script assigns lineage identity to ICM cells for the Spry4 dataset
## using a k-means clustering strategy 
## analogous to that we published in Saiz et al (2016) Nat Comms

## Function to do an exploratory plot of GATA6 vs NANOG values for ICM cells
spread.icm <- function(dataset) {
        ol <- ggplot(dataset, 
                     aes(x = CH5.ebLogCor, y = CH3.ebLogCor))
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

bb <- spry %>% filter(Treatment == 'Littermate', 
                      ## Exclude litters not stained for GATA6 and NANOG
                      !Litter %in% c('AG', 'V'), 
                      ## Do not subset for wild type embryos
                      ## as data shows a similar distribution for wt and hets
                      ## more data will yield denser clusters
                      #Genotype ==  'wt', 
                      Stage %in% c('32_64', '90_120', '120_150'),  
                      TE_ICM != 'TE') %>% 
        select(Cellcount, CH5.ebLogCor, CH3.ebLogCor)
## Scatter plot of selected data
spread.icm(bb)

set.seed(20170401)

## Calculate three centers for bb (DP, PrE and EPI)
pp <- kmeans(bb[, 2:3], 3)
## Extract the minimum center values for each channel 
## and append to the centers list as the DN center
dn <- c(pp$centers[, 1][which.min(pp$centers[, 1])], 
        pp$centers[, 2][which.min(pp$centers[, 2])])
centers <- rbind(pp$centers, dn, deparse.level = 0)

## Create conditions to test for ICM (is.icm) and morulas (is.morula)
is.icm <- spry$Stage != '<32' & spry$TE_ICM == 'ICM' & 
        !spry$Litter %in% c('V', 'AG')
is.te <- spry$Stage != '<32' & spry$TE_ICM == 'TE'
is.morula <- spry$Stage == '<32'

## Make matrix to hold sum of squares (rows = icm rows, 4 columns)
ssq <- matrix(0, length(spry$TE_ICM[is.icm]), 4)
## and populate with min sum of squares for each cell
## (min distance to each of the centers for each cell's GATA6 vs NANOG value)
for(i in 1:4) {
        ssq[,i] <- (spry$CH5.ebLogCor[is.icm] - centers[i,1])^2 + 
                (spry$CH3.ebLogCor[is.icm] - centers[i,2])^2
}
## calculate what center each cell is closest to 
## (which sum of squares is smallest)
min.ssq <- apply(ssq, 1, which.min)

## Create new variable in spry to hold the lineage identity
## assigned using k-means clustering (Identity.km)
spry$Identity.km <- rep(NA, nrow(spry))
## Morulas and TE cells remain unchaged
spry$Identity.km[is.morula] <- 'morula'
spry$Identity.km[is.te] <- 'TE'
## Assign identity to ICM cells based on min.ssq values
spry$Identity.km[is.icm] <- c('PRE', 'DP', 'EPI', 'DN')[min.ssq]
spry$Identity.km <- factor(spry$Identity.km, 
                           levels = c('DN', 'EPI', 'DP', 'PRE', 
                                      'TE', 'morula'))
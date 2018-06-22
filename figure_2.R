## Load standard plotting aesthetics
source('plotting-aes.R')
library('colorspace')

## Extract Nestor's data from main table
spry <- subset(spry, Experimenter == 'NS')

# Make a variable to score inside vs outside ('in' and 'ICM' vs 'out' and 'TE')
# solely for plotting purposes
spry$out_in <- ifelse(spry$TE_ICM %in% c('TE', 'out'), 'out', 'in')
spry$out_in <- factor(spry$out_in, levels = c('out', 'in'))

################################################################################
## Supplementary Figure 2a - Validation of ICM-specific expression #############
################################################################################

## Plot Venus (no.ab) or anti-GFP::AF488 levels 
## in ICM vs TE cells of het embryos, per stage, 
## compared to levels in wt and negative controls (secondary Ab only)
figs2a <- ggplot(data = spry %>% 
                         # Select wt and heterozygous littermates
                         filter(Treatment == 'Littermate',
                                Litter != 'Y', 
                                Genotype1 == 'het', 
                                Stage != '8_16'), 
                 aes(x = Stage, y = green.ebLogCor))
figs2a <- figs2a + geom_jitter(data = spry %>% 
                                       filter(Treatment == 'Littermate', 
                                              Litter != 'Y', 
                                              Genotype1 == 'wt', 
                                              Stage != '8_16'), 
                               color = I('gray65'), alpha = I(0.15), 
                               position = position_jitter(width = 0.35))
figs2a <- figs2a + geom_jitter(aes(color = out_in), 
                               size = 1.2, alpha = I(0.2), 
                               position = 
                                       position_jitterdodge(dodge.width = 0.9, 
                                                            jitter.width = 0.5))
figs2a <- figs2a + geom_boxplot(aes(fill = out_in), color = I('black'), 
                                position = position_dodge(width = 0.9), 
                                outlier.shape = 1, alpha = I(0))
figs2a <- figs2a + geom_jitter(data = spry %>%
                                       filter(Treatment == 'neg.control'), 
                               color = I('black'), alpha = I(0.15))
figs2a <- figs2a + scale_color_manual(values = idcols)
figs2a <- figs2a + scale_fill_manual(values = idcols)
figs2a <- figs2a + looks + facet_grid(green.ab2 ~ .) + coord_fixed(0.75)
figs2a <- figs2a + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(figs2a)

################################################################################
## Supplementary Figure 2b - Does the H2B-Venus casette affect embryo growth? ##
################################################################################

## Plot, for each stage and genotype, 
## the ratio of each embryo's cell count
## against the average cell count of its litter's heterozygotes
## (if there isn't any difference in size between genotypes
## this ratio should stay around 1 throughout)
rr <- subset(spry, Treatment == 'Littermate' &  
             Genotype1 != 'unknown' & 
             Cellcount > 15)
rr <- merge(rr, het.meancount)
rr$countratio <- rr$Cellcount / rr$het.meancount
rr <- rr %>% 
        group_by(Litter, Embryo_ID, Genotype1, countratio, Stage) %>% 
        summarize()

figs2b <- ggplot(data = rr, 
               aes(x = Stage, y = countratio)) 
figs2b <- figs2b + geom_boxplot(aes(fill = Genotype1), color = I('black'), 
                                outlier.shape = 1) 
figs2b <- figs2b + geom_jitter(aes(shape = Genotype1), 
                           color = I('black'), size = I(1.2),  
                           position = position_jitterdodge(dodge.width = 0.9, 
                                                           jitter.width = 0.5)) 
figs2b <- figs2b + ylim(0, 2) + scale_fill_manual(values = gencols)
figs2b <- figs2b + looks + coord_fixed(3)
print(figs2b)

## Plot total number of cells per litter, in ascending size
## comparing wild type (wt), Spry4:H2B-Venus/+ (het) 
## and Spry4:H2B-Venus/H2B-Venus embryos, for each litter
# figs2b <- ggplot(data = spry %>% 
#                          # Select embryos fixed upon collection (Littermates)
#                          # of known genotype
#                          filter(Treatment == 'Littermate', 
#                                 Genotype1 != 'unknown') %>% 
#                          group_by(Embryo_ID, 
#                                   Cellcount, 
#                                   litter.mean, 
#                                   Genotype1, 
#                                   Litter) %>% 
#                          summarize(), 
#                  aes(x = reorder(Litter, litter.mean), y = Cellcount))
# figs2b <- figs2b + geom_boxplot(aes(fill = Genotype1), color = I('black'), 
#                                 outlier.shape = 1)
# # figs2b <- figs2b + geom_violin(aes(fill = Genotype), color = I('black'),
# #                          scale = 'width')
# figs2b <- figs2b + geom_jitter(aes(shape = Genotype1), color = I('black'), 
#                                size = 2)
# figs2b <- figs2b + scale_fill_manual(values = gencols) 
# figs2b <- figs2b + labs(x = 'Litter', y = 'Total cell number', 
#                         fill = 'Genotype')
# figs2b <- figs2b + looks + coord_fixed(0.075)
# print(figs2b)

################################################################################
## Figure 2d - How do reporter levels (Venus) compare between genotypes? #######
################################################################################

# Make freq density plot for unstained embryos (H2B-Venus), 
# per stage, for ICM cells of each genotype and all TE cells
fig2d <- ggplot(data = spry %>% filter(Treatment == 'Littermate', 
                                       Genotype1 != 'unknown', 
                                       out_in == 'in', 
                                       green.marker == 'Venus', 
                                       Stage != '8_16'), 
                aes(x = green.ebLogCor))
fig2d <- fig2d + geom_density(data = spry %>% filter(Treatment == 'Littermate', 
                                                     Genotype1 != 'unknown', 
                                                     out_in == 'out', 
                                                     green.marker == 'Venus', 
                                                     Stage != '8_16'), 
                              aes(x = green.ebLogCor), fill = I('gray65'), 
                              color = I('black'), alpha = 0.75)
fig2d <- fig2d + geom_density(data = spry %>% filter(Treatment == 'Littermate', 
                                                     Genotype1 != 'unknown',
                                                     out_in == 'in', 
                                                     green.marker == 'Venus', 
                                                     Stage != '8_16'), 
                              aes(fill = Genotype1), color = I('black'), 
                              alpha = 0.75, trim = T)
fig2d <- fig2d + coord_fixed() + looks + scale_fill_manual(values = gencols)
fig2d <- fig2d + facet_grid(Stage ~ ., scales = 'free_y')
fig2d <- fig2d + theme(aspect.ratio = 0.25)
fig2d <- fig2d + labs(x = 'log[H2B-Venus]', y = 'Frequency', fill = 'Genotype')
print(fig2d)

################################################################################
## Figure 2e - How do reporter levels (anti-GFP) compare between genotypes? ####
################################################################################

# Calculate summary statistics of negative controls (no primary ab)
neg.con <- summary(subset(spry, Treatment == 'neg.control')$green.ebLogCor)
class(neg.con) <- 'numeric'

# Make freq density plot for embryos stained with anti-GFP::AF488, 
# per stage, for ICM cells of each genotype and all TE cells
fig2e <- ggplot(data = spry %>% filter(Treatment == 'Littermate', 
                                       Genotype1 != 'unknown', 
                                       out_in == 'in', 
                                       Litter != 'Y', 
                                       green.marker == 'GFP.ck', 
                                       Stage != '8_16'), 
                aes(x = green.ebLogCor))
fig2e <- fig2e + geom_density(data = spry %>% filter(Treatment == 'Littermate', 
                                                     Genotype1 != 'unknown', 
                                                     out_in == 'out', 
                                                     Litter != 'Y', 
                                                     green.marker == 'GFP.ck', 
                                                     Stage != '8_16'), 
                              aes(x = green.ebLogCor), fill = I('gray65'), 
                              color = I('black'), alpha = 0.75)
fig2e <- fig2e + geom_density(data = spry %>% filter(Treatment == 'Littermate', 
                                                      Genotype1 != 'unknown',
                                                     out_in == 'in', 
                                                     Litter != 'Y', 
                                                     green.marker == 'GFP.ck', 
                                                     Stage != '8_16'), 
                              aes(fill = Genotype1), color = I('black'), 
                              alpha = 0.75, trim = T)
fig2e <- fig2e + coord_fixed() + looks + scale_fill_manual(values = gencols)
fig2e <- fig2e + geom_vline(aes(xintercept = neg.con[3]), size = 0.75)
fig2e <- fig2e + geom_vline(aes(xintercept = neg.con[1]), linetype = 2)
fig2e <- fig2e + geom_vline(aes(xintercept = neg.con[6]), linetype = 2)
fig2e <- fig2e + facet_grid(Stage ~ ., scales = 'free_y') 
fig2e <- fig2e + theme(aspect.ratio = 0.25)
fig2e <- fig2e + labs(x = 'log[anti-GFP::AF488]', y = 'Frequency', 
                      fill = 'Genotype')
print(fig2e)

################################################################################
## Supplementary figure 2c - Centres for lineage assignment based on k-means ###
################################################################################

## Plot log[GATA6] vs log[NANOG] (corrected) for all ICM cells in littermates 
## of all genotypes, color coding for cell count (age) 
## and overlaying the cluster centers calculated in identify_spry.R
figs2c <- ggplot(data = spry %>% filter(Treatment == 'Littermate', 
                                        farred.marker == 'NANOG.rb', 
                                        red.marker == 'GATA6.gt',  
                                        Genotype1 != 'unknown', 
                                        TE_ICM %in% c('ICM', 'in'), 
                                        Stage != '<32'), 
                 aes(x = red.ebLogCor, y = farred.ebLogCor))
figs2c <- figs2c + geom_jitter(aes(color = Cellcount), size = 1.5)
figs2c <- figs2c + annotate('text', x = centers.ns[, 1], 
                            y = centers.ns[, 2], 
                            label =  c('DP', 'EPI', 'PrE', 'DN'), 
                            size = 10)
figs2c <- figs2c + looks + coord_fixed() + xlim(3, 8.25) + ylim(3, 8.25)
figs2c <- figs2c + labs(x = 'log[GATA6]', y = 'log[NANOG]', 
                        color = 'Cell count')
figs2c <- figs2c + scale_color_distiller(direction = 1, 
                                         palette = 'Blues')
print(figs2c)

################################################################################
## Supplementary figure 2d - Lineage assignment output per stage and genotype ##
################################################################################

## Plot log[GATA6] vs log[NANOG] (corrected) for all ICM cells in littermates
## separating by genotype and stage and color coding for identity
figs2d <- ggplot(data = spry %>% filter(Treatment == 'Littermate', 
                                        farred.marker == 'NANOG.rb', 
                                        red.marker == 'GATA6.gt',  
                                        Genotype1 != 'unknown',
                                        TE_ICM %in% c('ICM', 'in'), 
                                        Stage != '8_16'), 
                 aes(x = red.ebLogCor, y = farred.ebLogCor))
figs2d <- figs2d + geom_jitter(aes(color = Identity.km), size = 1.2)
figs2d <- figs2d + looks + coord_fixed() + facet_grid(Genotype1 ~ Stage)
figs2d <- figs2d + scale_color_manual(values = idcols) + 
        xlim(3, 8.25) + ylim(3, 8.25)
figs2d <- figs2d + labs(y = 'log[NANOG]', x = 'log[GATA6]', color = 'Identity')
print(figs2d)

################################################################################
## Figure 2f - How do reporter levels compare between ICM lineages #############
## in Spry4 hets? ##############################################################
################################################################################

## Plot reporter levels in each ICM lineage for each stage (> 15 cells)
## for ICM cells in heterozygous embryos only, both stainings, as boxplots
fig2f <- ggplot(data = spry %>% 
                           filter(Treatment == 'Littermate', 
                                  out_in == 'in', 
                                  Genotype1 != 'unknown', 
                                  Litter != 'Y', 
                                  farred.marker == 'NANOG.rb', 
                                  red.marker == 'GATA6.gt', 
                                  Stage != '8_16'), 
                   aes(x = Stage, y = green.ebLogCor))
fig2f <- fig2f + geom_jitter(data = spry %>% 
                                     filter(Treatment == 'Littermate', 
                                            TE_ICM != 'TE', 
                                            Genotype1 == 'wt', 
                                            Litter != 'Y', 
                                            farred.marker == 'NANOG.rb', 
                                            red.marker == 'GATA6.gt', 
                                            Stage != '8_16'),
                             color = I('gray65'), size = 1.2, 
                             alpha = 0.2, position = 
                                     position_jitter(width = 0.35))
fig2f <- fig2f + geom_boxplot(data = spry %>% 
                                            filter(Treatment == 'Littermate', 
                                                   out_in == 'in', 
                                                   Genotype1 == 'het', 
                                                   Litter != 'Y', 
                                                   farred.marker == 'NANOG.rb', 
                                                   red.marker == 'GATA6.gt', 
                                                   Stage != '8_16'), 
                                    aes(fill = Identity.km), 
                                    color = I('black'), outlier.shape = 1)
fig2f <- fig2f + looks + coord_fixed(1) + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
fig2f <- fig2f + facet_grid(green.ab2 ~ .) + labs(fill = 'Identity')
fig2f <- fig2f + scale_fill_manual(values = idcols)
print(fig2f)

################################################################################
## Supplementary figure 2e - Does the H2B-Venus casette knock-in affect ########
## the lineage composition of transgenic embryos? ##############################
################################################################################

## Stacked bar plot of TE vs ICM (broken down by lineage) composition per stage,
## grouped by genotype
figs2e <- ggplot(data = spry %>% filter(Cellcount > 31, 
                                        Genotype1 != 'unknown', 
                                        Treatment == 'Littermate', 
                                        farred.marker == 'NANOG.rb', 
                                        red.marker == 'GATA6.gt'), 
                 aes(x = Stage, fill = Identity.km))
figs2e <- figs2e + geom_bar(position = 'fill')
figs2e <- figs2e + looks + coord_fixed(5) + facet_grid( ~ Genotype1) + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs2e <- figs2e + scale_fill_manual(values = idcols) + labs(y = '% of total')
print(figs2e)

################################################################################
## Supplementary figure 2f - Does the H2B-Venus casette knock-in affect ########
## the lineage composition of the ICM of transgenic embryos? ###################
################################################################################

## Stacked bar plot of ICM composition per stage, grouped by genotype
figs2f <- ggplot(data = spry %>% filter(Treatment == 'Littermate', 
                                        Genotype1 != 'unknown', 
                                        TE_ICM == 'ICM', 
                                        Cellcount > 31, 
                                        farred.marker == 'NANOG.rb', 
                                        red.marker == 'GATA6.gt'), 
                 aes(x = Stage, fill = Identity.km))
figs2f <- figs2f + geom_bar(position = 'fill')
figs2f <- figs2f + looks + coord_fixed(5) + facet_grid( ~ Genotype1) + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs2f <- figs2f + scale_fill_manual(values = idcols) + labs(y = '% of ICM')
print(figs2f)

################################################################################
## Supplementary figure 2g - How do reporter levels compare between ############
## ICM lineages in Spry4 homozygous ko? ########################################
################################################################################

## Plot reporter levels in each ICM lineage for each stage (> 15 cells)
## for ICM cells in homozygous embryos only, both stainings, as boxplots
figs2g <- ggplot(data = spry %>% 
                        filter(Treatment == 'Littermate', 
                               out_in == 'in', 
                               Genotype1 != 'unknown', 
                               Litter != 'Y', 
                               farred.marker == 'NANOG.rb', 
                               red.marker == 'GATA6.gt', 
                               Stage != '8_16'), 
                aes(x = Stage, y = green.ebLogCor))
figs2g <- figs2g + geom_jitter(data = spry %>% 
                                     filter(Treatment == 'Littermate', 
                                            TE_ICM != 'TE', 
                                            Genotype1 == 'wt', 
                                            Litter != 'Y', 
                                            farred.marker == 'NANOG.rb', 
                                            red.marker == 'GATA6.gt', 
                                            Stage != '8_16'),
                             color = I('gray65'), size = 1.2, 
                             alpha = 0.2, position = 
                                     position_jitter(width = 0.35))
figs2g <- figs2g + geom_boxplot(data = spry %>% 
                                      filter(Treatment == 'Littermate', 
                                             out_in == 'in', 
                                             Genotype1 == 'homo', 
                                             Litter != 'Y', 
                                             farred.marker == 'NANOG.rb', 
                                             red.marker == 'GATA6.gt', 
                                             Stage != '8_16'), 
                              aes(fill = Identity.km), 
                              color = I('black'), outlier.shape = 1)
figs2g <- figs2g + looks + coord_fixed(1) + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs2g <- figs2g + facet_grid(green.ab2 ~ .) + labs(fill = 'Identity')
figs2g <- figs2g + scale_fill_manual(values = idcols)
print(figs2g)

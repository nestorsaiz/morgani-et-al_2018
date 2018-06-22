## Load standard plotting aesthetics
source('plotting-aes.R')
library('colorspace')

## Extract Nestor's data from main table
spry <- subset(spry, Experimenter == 'NS')

################################################################################
## Supplementary Figure 3a - How does culture/treatment affect embryo growth? ##
################################################################################

## Plot total cell numbers for each experimental group
## (Control vs each treatment) as boxplots
figs3a <- ggplot(data = spry %>% 
                         filter(!Treatment %in% c('Littermate', 
                                                  'neg.control'), 
                                Genotype1 != 'unknown', 
                                Tt_stage == '8cell') %>% 
                         group_by(Embryo_ID, Cellcount, 
                                  Treatment, Litter,
                                  Gene1, Genotype1, 
                                  Gene2, Genotype2, 
                                  Experiment, green.marker) %>% 
                         summarize(), 
                 aes(x = Treatment, y = Cellcount))
figs3a <- figs3a + geom_boxplot(color = I('black'), fill = I('gray75'), 
                                outlier.shape = 1, outlier.size = 2)
figs3a <- figs3a + geom_jitter(aes(color = Genotype1), size = 1.5, 
                               position = 
                                       position_jitterdodge(dodge.width = 0.75, 
                                                            jitter.width = 0.5))
figs3a <- figs3a + looks + coord_fixed(0.075)
figs3a <- figs3a + scale_color_manual(values = gencols)
figs3a <- figs3a + theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs3a <- figs3a + labs(color = 'Genotype', y = 'Total cell number')
print(figs3a)

################################################################################
## Supplementary Figure 3b - Effect of culture/treatment on ICM gene expression#
################################################################################

## Scatter plots of log[GATA6] vs log[NANOG] (corrected) for all ICM cells
## for each treatment and genotype, color coded for identity
figs3b <- ggplot(data = spry %>% filter(!Treatment %in% c('Littermate', 
                                                          'neg.control'), 
                                        Genotype1 != 'unknown', 
                                        TE_ICM == 'ICM', 
                                        Tt_stage == '8cell'), 
                 aes(x = red.ebLogCor, y = farred.ebLogCor))
figs3b <- figs3b + geom_jitter(aes(color = Identity.km), size = 1.5)
figs3b <- figs3b + looks + coord_fixed() + facet_grid(Genotype1 ~ Treatment)
figs3b <- figs3b + scale_color_manual(values = idcols) 
figs3b <- figs3b + xlim(3.5, 8.5) + ylim(3.5, 8.5)
figs3b <- figs3b + labs(y = 'log[NANOG]', x = 'log[GATA6]', color = 'Identity')
print(figs3b)

################################################################################
## Supplementary Figure 3c - Effect of culture/treatment on ICM composition ####
################################################################################

## Plot ICM composition as % of the ICM, 
## for each treatment, grouped by genotype
figs3c <- ggplot(data = spry %>% filter(!Treatment %in% c('Littermate', 
                                                          'neg.control'), 
                                        Genotype1 != 'unknown', 
                                        TE_ICM == 'ICM', 
                                        Tt_stage == '8cell'), 
                 aes(x = Treatment, fill = Identity.km))
figs3c <- figs3c + geom_bar(position = 'fill')
figs3c <- figs3c + facet_grid( ~ Genotype1)
figs3c <- figs3c + looks + coord_fixed(8)
figs3c <- figs3c + theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs3c <- figs3c + scale_fill_manual(values = idcols) + 
        labs(y = '% of ICM', fill = 'Identity')
print(figs3c)

################################################################################
## Figure 3b - How does the treatment affect reporter expression? ##############
################################################################################

## Plot Venus levels (anti-GFP::AF488) in the ICM vs TE of heterozygous embryos
## for each treatment condition
## show the level in the corresponding wild type embryos for reference
fig3c <- ggplot(data = spry %>% filter(!Treatment %in% c('Littermate', 
                                                         'neg.control'), 
                                       Genotype1 == 'het', 
                                       TE_ICM == 'ICM', 
                                       Tt_stage == '8cell'), 
                aes(x = Treatment, y = green.ebLogCor))
fig3c <- fig3c + geom_jitter(data = spry %>% 
                                     filter(!Treatment %in% c('Littermate', 
                                                              'neg.control'),
                                            Genotype1 == 'wt', 
                                            Tt_stage == '8cell', 
                                            TE_ICM == 'ICM'), 
                             aes(color = Genotype1), size = 1.5, alpha = I(0.85), 
                             position = position_jitter(width = 0.35))
fig3c <- fig3c + geom_jitter(aes(color = Genotype1), size = 1.5, alpha = I(0.5), 
                             position = 
                                     position_jitterdodge(jitter.width = 1.8))
fig3c <- fig3c + geom_boxplot(aes(fill = Genotype1), color = I('black'), 
                              outlier.shape = 1, outlier.size = 2, alpha = I(0),
                              position = position_dodge(width = 0.9))
fig3c <- fig3c + geom_hline(aes(yintercept = neg.con[3]), size = 0.75)
fig3c <- fig3c + geom_hline(aes(yintercept = neg.con[1]), linetype = 2)
fig3c <- fig3c + geom_hline(aes(yintercept = neg.con[6]), linetype = 2)
fig3c <- fig3c + looks + coord_fixed(1.25) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
fig3c <- fig3c + scale_fill_manual(values = gencols)
fig3c <- fig3c + scale_color_manual(values = gencols)
fig3c <- fig3c + labs(y = 'log[anti-GFP::AF488]', fill = 'Genotype', 
                      color = 'Genotype')
print(fig3c)

################################################################################
## Figure 3e - Timeline of changes in reporter expression after FGF inhibition #
################################################################################

fig3e <- ggplot(data = spry.mov, 
                aes(x = h4, y = meanorm))
fig3e <- fig3e + geom_boxplot(aes(fill = Genotype1), color = 'black',  
                              outlier.shape = 1, outlier.size = 2)
# fig3e <- fig3e + geom_jitter(aes(color = Slice), size = 0.75, 
#                              position = position_jitter(width = 0.3), 
#                              alpha = 0.5)
#fig3e <- fig3e + scale_color_distiller(direction = 1, palette = 'Blues')
fig3e <- fig3e + scale_fill_manual(values = gencols)
fig3e <- fig3e + facet_wrap( ~ Treatment) + looks + coord_fixed(4)
fig3e <- fig3e + theme(axis.text.x = element_text(angle = 30, hjust = 1))
fig3e <- fig3e + labs(x = 'Hour bins', 
                      y = 'Normalized average intensity per z-slice, per hour')
print(fig3e)

################################################################################
## Figure 3h - Reporter expression levels in late blastocysts after ############
## FGF inhibition from mid blastocyst ##########################################
################################################################################

## Extract embryos cultured from blastocyst stage
spy.late <- subset(spry, Treatment %in% c('Control', 'AZD_1', 'PD03_1') & 
                           Tt_length != '48h')
spy.late$Treatment <- factor(spy.late$Treatment, 
                             levels = c('Littermate', 'Control', 
                                        'FGF4_1000', 'AZD_1', 
                                        'PD03_1', 'neg.control'))

fig3g <- ggplot(data = spy.late %>% filter(Litter != 'Y', 
                                        Genotype1 == 'het'), 
             aes(x = Treatment, y = green.ebLogCor))
fig3g <- fig3g + geom_jitter(aes(color = TE_ICM), size = I(1.5), alpha = I(0.75), 
                       position = 
                               position_jitterdodge(dodge.width = 0.9, 
                                                    jitter.width = 0.8))
fig3g <- fig3g + geom_boxplot(aes(fill = TE_ICM), color = I('black'), 
                        alpha = I(0), outlier.shape = 1, 
                        outlier.size = 2, 
                        position = position_dodge(width = 0.9))
fig3g <- fig3g + looks + scale_color_manual(values = idcols) + coord_fixed()
fig3g <- fig3g + labs(y = 'log[anti-GFP::AF488]', color = 'Identity', 
                      fill = 'Identity')
fig3g <- fig3g + ylim(4, 9)
print(fig3g)

################################################################################
## Supplementary Figure 3d - How does late treatment affect embryo growth? #####
################################################################################

figs3d <- ggplot(data = spy.late %>% 
                         group_by(Embryo_ID, Cellcount, Treatment, 
                                  Tt_stage, Tt_length, Litter, 
                                  icm.count, litter.mean, Stage, 
                                  Gene1, Genotype1) %>% 
                         summarize(), 
                 aes(x = Treatment, y = Cellcount))
figs3d <- figs3d + geom_boxplot(color = I('black'), fill = I('gray75'), 
                                outlier.shape = 1, outlier.size = 2)
figs3d <- figs3d + geom_jitter(aes(color = Genotype1), size = 1.5, 
                               position = 
                                       position_jitterdodge(dodge.width = 0.75, 
                                                            jitter.width = 0.5))
figs3d <- figs3d + scale_color_manual(values = gencols)
figs3d <- figs3d + looks + coord_fixed(0.1) + ylim(40, 120)
figs3d <- figs3d + theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs3d <- figs3d + labs(color = 'Genotype', y = 'Total cell number')
print(figs3d)

################################################################################
## Supplementary Figure 3e - Effect of late treatment on ICM gene expression ###
################################################################################

figs3e <- ggplot(data = spy.late %>% 
                         filter(TE_ICM == 'ICM'), 
                 aes(x = red.ebLogCor, y = farred.ebLogCor))
figs3e <- figs3e + geom_jitter(aes(color = Identity.km), size = 1.5)
figs3e <- figs3e + looks + coord_fixed() + facet_grid(Genotype1 ~ Treatment)
figs3e <- figs3e + scale_color_manual(values = idcols) 
figs3e <- figs3e + labs(y = 'log[NANOG]', x = 'log[GATA6]', color = 'Identity')
print(figs3e)

################################################################################
## Supplementary Figure 3f - Effect of late treatment on ICM composition #######
################################################################################

## Plot ICM composition as % of the ICM, 
## for each treatment, grouped by genotype
figs3f <- ggplot(data = spy.late %>% 
                         filter(TE_ICM == 'ICM', 
                                Genotype1 == 'het'), 
                 aes(x = Treatment, fill = Identity.km))
figs3f <- figs3f + geom_bar(position = 'fill')
figs3f <- figs3f + looks + coord_fixed(6)
figs3f <- figs3f + theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs3f <- figs3f + scale_fill_manual(values = idcols) + 
        labs(y = '% of ICM', fill = 'Identity')
print(figs3f)

################################################################################
## Supplementary Figure 3g - Reporter expression in each ICM lineage ###########
## of late blastocysts after FGF inhibition from mid blastocyst ################
################################################################################

figs3g <- ggplot(data = spy.late %>% filter(TE_ICM == 'ICM', 
                                         Litter != 'Y', 
                                         Genotype1 == 'het'), 
              aes(x = Treatment, y = green.ebLogCor))
figs3g <- figs3g + geom_jitter(data = spy.late %>% filter(TE_ICM == 'TE', 
                                                    Litter != 'Y', 
                                                    Genotype1 == 'het'), 
                         aes(color = TE_ICM), 
                         alpha = I(0.5), 
                         position = position_jitter(width = 0.35))
figs3g <- figs3g + geom_boxplot(aes(fill = Identity.km), color = I('black'), 
                          outlier.shape = 1, outlier.size = 2, 
                          position = position_dodge(width = 0.9))
figs3g <- figs3g + looks + coord_fixed()
figs3g <- figs3g + scale_color_manual(values = idcols) 
figs3g <- figs3g + scale_fill_manual(values = idcols)
figs3g <- figs3g + labs(y = 'log[anti-GFP::AF488]', color = 'Identity', 
                      fill = 'Identity')
print(figs3g)

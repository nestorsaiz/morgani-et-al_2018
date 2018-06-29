source('plotting-aes.R')
spry_vg <- dplyr::filter(spry, Experimenter == 'VG')
spry_vg$Identity.km <- factor(spry_vg$Identity.km, levels = c('TE', 'PRE',
                                                              'DP', 'EPI',
                                                              'DN', 'morula'))

spry_vg$Genotype2 <- as.character(spry_vg$Genotype2)
spry_vg$Genotype22 <- interaction(spry_vg$Gene2, spry_vg$Genotype2)
is.wt <- spry_vg$Genotype22 %in% c('Fgf4.wt', 'Fgfr1.wt', 'Fgfr2.wt')
spry_vg$Genotype22 <- as.character(spry_vg$Genotype22)
spry_vg$Genotype22[is.wt] <- 'wt'
spry_vg$Genotype22 <- factor(spry_vg$Genotype22, levels = c('wt', 'Fgf4.ko',
                                                            'Fgf4.het',
                                                            'Fgfr1.ko',
                                                            'Fgfr1.het',
                                                            'Fgfr2.ko',
                                                            'Fgfr2.het'))
spry_vg_sub <- dplyr::filter(spry_vg, Genotype2 %in% c('wt', 'ko'),
                               Stage != '120_150')


## Calculate the average level for each fluorescence channel
## for each embryo and lineage
meantensity <- spry %>%
  group_by(Experiment, Litter, Embryo_ID, Identity.km, TE_ICM, 
           Cellcount, Stage, Exp_date, Img_date, Genotype) %>%
  summarize(Count = n(), 
            CH1.mean = mean(CH1.ebLogCor),
            CH2.mean = mean(CH2.ebLogCor), 
            CH3.mean = mean(CH3.ebLogCor),
            CH4.mean = mean(CH4.ebLogCor),
            CH5.mean = mean(CH5.ebLogCor))

## Calculate the number of cells per lineage for each embryo
spry.lincounts <- spry %>% 
  group_by(Experiment, Litter, 
           Embryo_ID, TE_ICM, 
           Cellcount, Stage, 
           Exp_date, Img_date, 
           Genotype, Identity.km) %>%
  summarize(count = n())

################################################################################
################################################################################

## Main figure plots
#Figure 4
fig4 <- ggplot(data = subset(spry_vg, Identity.km == 'EPI' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150' | 
                             Identity.km == 'PRE' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150' | 
                             Identity.km == 'DP' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150'), 
             aes(x = Genotype22, y = green.ebLogCor))
fig4 <- fig4 + geom_jitter(aes(color = Identity.km),
                       size = 1.2, #alpha = I(0.5), 
                       position = 
                         position_jitterdodge(jitter.width = 0.65))
fig4 <- fig4 + geom_boxplot(aes(fill = Identity.km), color = I('black'),
                        alpha = 0.0, 
                        outlier.shape = 1)
fig4 <- fig4 + scale_color_manual(values = idcols)
fig4 <- fig4 + scale_fill_manual(values = idcols)
fig4 <- fig4 + looks + facet_grid( ~Stage) + coord_fixed(0.75)
fig4 <- fig4 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(fig4)

#fig4 <- fig4 + geom_jitter(data = subset(spry_vg, Identity.km == 'TE' &
                                       #Genotype22 != 'Fgf4.het' &
                                       #Genotype22 != 'Fgfr1.het' &
                                       #Genotype22 != 'Fgfr2.het' &
                                       #Stage != '120_150'),
                       #aes(fill = Identity.km), color = 'green',
                       #size = 0.5, alpha = I(0.2))

## Supplementary figure plots
# Figure S4
figs4 <- ggplot(data = subset(spry_vg, TE_ICM == 'TE' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150' | 
                             TE_ICM == 'ICM' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150'), 
             aes(x = Genotype22, y = green.ebLogCor))
figs4 <- figs4 + geom_jitter(aes(color = TE_ICM),
                       size = 1.2, #alpha = I(0.2), 
                       position = 
                         position_jitterdodge(jitter.width = 0.65))
figs4 <- figs4 + geom_boxplot(aes(fill = TE_ICM), color = I('black'),
                        alpha = 0.0, 
                        outlier.shape = 1)
figs4 <- figs4 + scale_color_manual(values = idcols)
figs4 <- figs4 + scale_fill_manual(values = idcols)
figs4 <- figs4 + looks + facet_grid( ~Stage) + coord_fixed(0.75)
figs4 <- figs4 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(figs4)


figs4a <- ggplot(data = subset(spry_vg, TE_ICM == 'ICM' &
                                 Genotype22 != 'Fgf4.het' &
                                 Genotype22 != 'Fgfr1.het' &
                                 Genotype22 != 'Fgfr2.het' &
                                 Stage != '120_150'),
                 aes(x = red.ebLogCor, y = farred.ebLogCor))
figs4a <- figs4a + geom_jitter(aes(color = Identity.km), size = 1.2)
figs4a <- figs4a + looks + coord_fixed() + facet_grid(Genotype22 ~ Stage)
figs4a <- figs4a + scale_color_manual(values = idcols) + 
  xlim(3, 9) + ylim(3, 9)
figs4a <- figs4a + labs(y = 'log[NANOG]', x = 'log[GATA6]', color = 'Identity')
print(figs4a)

figs4b <- ggplot(data = subset(spry_vg, TE_ICM == 'ICM' &
                                 Genotype22 != 'Fgf4.het' &
                                 Genotype22 != 'Fgfr1.het' &
                                 Genotype22 != 'Fgfr2.het' &
                                 Stage != '120_150'),
                 aes(x = Stage, fill = Identity.km))
figs4b <- figs4b + geom_bar(position = 'fill')
figs4b <- figs4b + looks + coord_fixed(5) + facet_grid( ~ Genotype22) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
figs4b <- figs4b + scale_fill_manual(values = idcols) + labs(y = '% of ICM')
print(figs4b)


################################################################################
################################################################################
## Miscellaneous plots
p3 <- ggplot(data = subset(spry_vg, TE_ICM == 'TE' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150'), 
             aes(x = Genotype22, y = (Cellcount - icm.count)))
p3 <- p3 + geom_jitter(aes(color = TE_ICM),
                       size = 1.2, #alpha = I(0.2), 
                       position = 
                         position_jitterdodge(dodge.width = 0.9, 
                                              jitter.width = 0.65))
p3 <- p3 + geom_boxplot(aes(fill = TE_ICM), color = I('black'),
                        alpha = 0.0,
                        position = position_dodge(width = 0.9), 
                        outlier.shape = 1)
p3 <- p3 + scale_color_manual(values = idcols)
p3 <- p3 + scale_fill_manual(values = idcols)
p3 <- p3 + looks + facet_grid( ~Stage)
p3 <- p3 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p3)

p4 <- ggplot(data = subset(spry_vg, TE_ICM == 'ICM' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150'), 
             aes(x = Genotype22, y = icm.count))
p4 <- p4 + geom_jitter(aes(color = TE_ICM),
                       size = 1.2, #alpha = I(0.2), 
                       position = 
                         position_jitterdodge(dodge.width = 0.9, 
                                              jitter.width = 0.65))
p4 <- p4 + geom_boxplot(aes(fill = TE_ICM), color = I('black'),
                        alpha = 0.0,
                        position = position_dodge(width = 0.9), 
                        outlier.shape = 1)
p4 <- p4 + scale_color_manual(values = idcols)
p4 <- p4 + scale_fill_manual(values = idcols)
p4 <- p4 + looks + facet_grid( ~Stage)
p4 <- p4 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p4)

p5 <- ggplot(data = subset(spry_vg, Identity.km == 'EPI' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150'| 
                             Identity.km == 'PRE' &
                             Genotype22 != 'Fgf4.het' &
                             Genotype22 != 'Fgfr1.het' &
                             Genotype22 != 'Fgfr2.het' &
                             Stage != '120_150'), 
             aes(x = Genotype22, y = green.ebLogCor))
p5 <- p5 + geom_jitter(aes(fill = Identity.km), color = I('black'), 
                       size = 1.2, alpha = I(0.2), 
                       position = 
                         position_jitterdodge(dodge.width = 0.9, 
                                              jitter.width = 0.5))
p5 <- p5 + geom_jitter(data = subset(spry_vg, Identity.km == 'TE' &
                                       Genotype22 != 'Fgf4.het' &
                                       Genotype22 != 'Fgfr1.het' &
                                       Genotype22 != 'Fgfr2.het' &
                                       Stage != '120_150'),
                       aes(fill = Identity.km), color = 'green',
                       size = 0.5, alpha = I(0.1))
p5 <- p5 + geom_boxplot(aes(fill = Identity.km), color = I('black'), 
                        position = position_dodge(width = 0.9), 
                        outlier.shape = 1)
p5 <- p5 + scale_color_manual(values = idcols)
p5 <- p5 + scale_fill_manual(values = idcols)
p5 <- p5 + looks + facet_grid( ~Stage) + coord_fixed(0.75)
p5 <- p5 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p5)

## Miscellaneous lines to test & extract data

# Extracts Embryo_ID of particular genotype and stage
unique(subset(spry, Gene2 == 'Fgf4' &
                Genotype2 == 'wt' & Stage == '90_120')$Embryo_ID)

#Plotting individual cells by expression levels of
#Farred/Red channels and color-coded by Green channel;
#faceted by stage & genotype

qplot(red.ebLogCor,  farred.ebLogCor,
      data = subset(spry_vg, TE_ICM == 'ICM' & Treatment == "Littermate" & 
                      !Litter %in% c('AG', 'V') & Genotype1 != 'unknown' & 
                      Genotype2  != 'het'), 
      color = green.ebLogCor) + looks + #scale_color_manual(values = idcols) + 
  facet_grid(Genotype2 + Gene2 ~ Stage) + coord_fixed()
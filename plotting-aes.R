library('ggplot2')

## Create vector to define lineage identity colors
idcols <- c('EPI' = 'red', 'PRE' = 'blue', 'DP' = 'purple', 
            'DN' = 'gray', 'TE' = 'green', 'ICM' = 'purple', 
            'morula' = 'violetred', 'ESC' = 'red4', 'out' = 'green', 
            'in' = 'purple')

## Create vector to define cell type color (donor vs host)
cellcols <- c('host' = 'gray', 'donor' = 'green')

## And genotype colors
gencols <- c('wt' = '#D6F9DD', 'het' = '#379A54', 
             'ko' = '#007300', 'unknown' = 'black')

## Make object containing aesthetics for the plots (font size, etc)
looks <- theme_bw() + theme(axis.text = element_text(size = 15, 
                                                     color = 'black'), 
                            axis.title = element_text(size = 20, 
                                                      color = 'black'), 
                            legend.text = element_text(size = 15, 
                                                       color = 'black'), 
                            legend.title = element_text(size = 20, 
                                                        color = 'black'), 
                            axis.ticks = element_line(color = 'black', 
                                                      size = 0.5), 
                            panel.grid = element_blank(), 
                            strip.background = element_blank(), 
                            strip.text = element_text(size = 17.5, 
                                                      color = 'black'), 
                            panel.border = element_rect(color = 'black', 
                                                        size = 1))
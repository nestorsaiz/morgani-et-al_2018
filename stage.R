stage <- function(dataset){
        ## Given a certain dataset, stage embryos according to their cell count
        ## 1. Define the intervals to stage in (cuts)
        cuts <- c(1, 8, 16, 32, 64, 90, 120, 150, Inf)
        ## 2. Label the stages (based on intervals) as desired
        stages <- c('<8', '8_16', '16_32', '32_64', '64_90', 
                    '90_120', '120_150', '>150')
        ## 3. Apply the cut function for the defined intervals and labels
        ## leaving intervals open on the right (excludes upper limit)
        dataset$Stage <- cut(dataset$Cellcount, breaks = cuts, 
                             labels = stages, right = F)
        # Convert 'Stage' into a factor with the levels ordered
        # in increasing number of cells
        dataset$Stage <- factor(dataset$Stage, levels = c('<8', '8_16', '16_32', 
                                                          '32_64', '64_90', 
                                                          '90_120', '120_150', 
                                                          '>150'))
        return(dataset)
}

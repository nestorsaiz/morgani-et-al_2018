## Create a function that calculates the mode for Channel 2 (Venus or anti-GFP)
## for ICM cells at each stage of Spry4:H2B-Venus/+ ('het') embryos 
## fixed on collection
moda <- function(dataset, stage, marker) {
        # Subset as described above
        x <- subset(dataset, Treatment == 'Littermate' & 
                            Genotype == 'het' & TE_ICM != 'TE' & 
                            !Litter %in% c('V', 'AG') & 
                            Stage == stage & venus.gfp == marker)$CH2.ebLogCor
        # Calculate the density distribution of Channel 2 values (x)
        d <- density(x)
        # Extract the mode and the maximum value
        c(d$x[which.max(d$y)], max(d$y))
}

## Create a matrix to be filled with the mode and max for each stage

# Extract the levels that Stage takes in the dataset
spry.stage <- levels(spry$Stage)
#spry.stage <- spry.stage[2:6]
# Create a matrix with columns for Stage, Mode and Max
modas <- matrix(c(1:(3*length(spry.stage))), nrow = length(spry.stage), 
                ncol = 3,
                dimnames = list(c(), c('Stage', 'Mode', 'Max')))
# Fill column 1 with the levels of Stage
modas[, 1] <- spry.stage
# For each stage in unstained embryos (CH2 = 'Venus'), calculate the mode and max
# and store it in the corresponding row in modas
for(s in spry.stage) {
        i <- which(spry.stage == s)
        modas[i, 2:3] <- moda(spry, stage = s, marker = 'Venus')[1:2]
}
# Duplicate modas to moda.venus
moda.venus <- as.data.frame(modas, stringsAsFactors = FALSE)
moda.venus$Mode <- as.numeric(moda.venus$Mode)
moda.venus$Max <- as.numeric(moda.venus$Max)

# Repeat the operation above for embryos stained with anti-GFP (CH2 = 'GFP.ck')
for(s in spry.stage) {
        i <- which(spry.stage == s)
        modas[i, 2:3] <- moda(spry, stage = s, marker = 'GFP.ck')[1:2]
}
# Duplicate modas to moda.gfp
moda.gfp <- as.data.frame(modas, stringsAsFactors = FALSE)
moda.gfp$Mode <- as.numeric(moda.gfp$Mode)
moda.gfp$Max <- as.numeric(moda.gfp$Max)
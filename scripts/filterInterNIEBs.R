### OUTPUT SPECIAL BED-LIKE FORMAT
library(readr)


args = commandArgs(trailingOnly=TRUE) # input, output
# args = c("/home/fabien/Documents/20200228_interSmallNFR.dat", "test_output.pouet")
data = read_tsv(args[1], col_names = FALSE)

# data = read_tsv("../Data/20200228_interSmallNFR.dat", col_names = FALSE)
# data = data[data$X4 - data$X3 > 1000,]
data$middle = as.integer((data$X3 + data$X4) / 2)
data$BL = as.integer((data$X2 + data$X3) /2)
data$ER = as.integer((data$X4 + data$X5) / 2)

Left = data.frame("chrom" = data$X1, "Begin" = data$BL, "Zero" = data$X3, "End" = data$middle, "Sense" = "L")
Right = data.frame("chrom" = data$X1, "Begin" = data$middle, "Zero" = data$X4, "End" = data$ER, "Sense" = "R")

output = rbind(Left, Right)
output = output[order(output$chrom, output$Begin),]

## put it in bed form
output$name = "."
output$score = "."

output = output[,c(1, 2, 4, 6, 7, 5, 3)]

write_tsv(output, args[2], col_names = FALSE)

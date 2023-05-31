### OUTPUT SPECIAL BED-LIKE FORMAT
# library(readr)

options(scipen=999)
args = commandArgs(trailingOnly=TRUE) # input, output
# args = c("/home/fabien/Documents/20200228_interSmallNFR.dat", "test_output.pouet")
data = read.delim(args[1], header = F)

if(length(args) > 2) {
    min_intra = as.numeric(args[3])
    min_extra = as.numeric(args[4])
    print(min_intra)
} else {
    min_extra = 0
    min_intra = 0
}

# data = read_tsv("../Data/20200228_interSmallNFR.dat", col_names = FALSE)
# data = data[data$X4 - data$X3 > 1000,]

Left = data.frame("chrom" = data$V1, "start" = ceiling((data$V2 + data$V3) / 2), "end" = floor((data$V3 + data$V4) / 2), "strand" = "L", "zero" = (data$V3))
Right = data.frame("chrom" = data$V1, "start" = ceiling((data$V3 + data$V4) / 2), "end" = floor((data$V4 + data$V5) / 2), "strand" = "R", "zero" = (data$V4))

# data$middle = round((data$V3 + data$V4) / 2)
# data$BL = round((data$V2 + data$V3) / 2)
# data$ER = round((data$V4 + data$V5) / 2)

# Left = data.frame("chrom" = data$X1, "start" = data$BL, "zero" = data$X3 - 1, "end" = data$middle, "strand" = "L") # correct for half open interval
# Right = data.frame("chrom" = data$X1, "start" = data$middle, "zero" = data$X4 - 1, "end" = data$ER, "strand" = "R")

Left = Left[(Left$zero - Left$start) > (min_intra * -1),]
Left = Left[(Left$end - Left$zero) > min_extra,]

Right = Right[(Right$zero - Right$start) > min_extra,]
Right = Right[(Right$end - Right$zero) > (min_intra * -1),]

output = rbind(Left, Right)
output = output[order(output$chrom, output$start),]

## put it in bed form
output$name = "."
output$score = "."

output = output[,c("chrom", "start", "end", "name", "score", "strand", "zero")]

write.table(output, file = args[2], col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

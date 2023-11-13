### OUTPUT SPECIAL BED-LIKE FORMAT

options(scipen = 999)
args = commandArgs(trailingOnly = TRUE) # input, output
# args = c("/home/fabien/Documents/20200228_interSmallNFR.dat", "test_output.pouet")
data = read.delim(args[1], header = FALSE)

if(length(args) > 2) {
    min_intra = as.numeric(args[3])
    min_extra = as.numeric(args[4])
} else {
    min_extra = 0
    min_intra = 0
}

if(length(args) > 4) {
    max_intra = as.numeric(args[5])
    max_extra = as.numeric(args[6])
} else {
    max_intra = Inf
    max_extra = Inf
}

print(min_extra)
print(max_extra)
print(min_intra)
print(max_intra)
# data = read_tsv("../Data/20200228_interSmallNFR.dat", col_names = FALSE)
# data = data[data$X4 - data$X3 > 1000,]

left = data.frame("chrom" = data$V1,
                  "start" = ceiling((data$V2 + data$V3) / 2),
                  "end" = floor((data$V3 + data$V4) / 2),
                  "strand" = "+", "zero" = (data$V3))
right = data.frame("chrom" = data$V1,
                   "start" = ceiling((data$V3 + data$V4) / 2),
                   "end" = floor((data$V4 + data$V5) / 2),
                   "strand" = "-", "zero" = (data$V4))
intra_right = right$end - right$zero
extra_right = right$zero - right$start
valid_right = intra_right > min_intra & intra_right < max_intra & extra_right > min_extra & extra_right < max_extra
intra_left = left$zero - left$start
extra_left = left$end - left$zero
valid_left = intra_left > min_intra & intra_left < max_intra & extra_left > min_extra & extra_left < max_extra

left = left[valid_left,]
right = right[valid_right,]

output = rbind(left, right)
output = output[order(output$chrom, output$start), ]

## put it in bed form
output$name = "."
output$score = "."

output = output[, c("chrom", "start", "end", "name", "score", "strand", "zero")]

write.table(output, file = args[2], col.names = FALSE,
            row.names = FALSE, quote = FALSE, sep = "\t")

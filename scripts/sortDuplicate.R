#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

print("Filtering overlaps...")
overlaps = read.delim(args[1], header = FALSE)
source = read.delim(args[2], header = FALSE)

not_overlaping = overlaps[overlaps$V7 < 2, ]
not_overlaping = not_overlaping$V4

write.table(source[source$V4 %in% not_overlaping, ],
            file = args[3], col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
write.table(source[!(source$V4 %in% not_overlaping), ],
            file = args[4], col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

count_overlaps = nrow(source) - nrow(source[source$V4 %in% not_overlaping, ])

print("Finished filtering.")
print(paste("Filtered", count_overlaps, "rows"))

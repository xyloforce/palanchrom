#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

print("Filtering overlaps...")
overlaps = read.table(args[1], sep = "\t", header = FALSE)
source = read.table(args[2], sep = "\t", header = FALSE)
not_overlaping = overlaps[overlaps$V7 < 2, "V4"]
write.table(source[source$V4 %in% not_overlaping, ], file = args[3], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(source[!(source$V4 %in% not_overlaping), ], file = args[4], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("Finished filtering.")
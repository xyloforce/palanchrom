#!/usr/bin/env Rscript
library("readr")

args = commandArgs(trailingOnly = TRUE)

print("Filtering overlaps...")
overlaps = read_tsv(args[1], col_names = FALSE, show_col_types = FALSE)
source = read_tsv(args[2], col_names = FALSE, show_col_types = FALSE)

not_overlaping = overlaps[overlaps$X7 < 2,]
not_overlaping = not_overlaping$X4

write_tsv(source[source$X4 %in% not_overlaping, ], file = args[3], col_names = FALSE)
write_tsv(source[!(source$X4 %in% not_overlaping), ], file = args[4], col_names = FALSE)

count_overlaps = nrow(source) - nrow(source[source$X4 %in% not_overlaping, ])

print("Finished filtering.")
print(paste("Filtered", count_overlaps, "rows"))

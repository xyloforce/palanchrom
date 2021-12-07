library(readr)


args = commandArgs(trailingOnly=TRUE) # input, output
data = read_tsv(args[1], col_names = TRUE, skip = 1)

data = data[!(toupper(data$ALT) == "N"),]
head(data)
write_tsv(data, args[2])

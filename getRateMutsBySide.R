library(readr)
library(ggplot2)
library(extrafont)

args = commandArgs(trailingOnly=TRUE)

source("~/setThemePoster.R")

reverseBase = function(x) {
	return(switch(x[["base"]], "A"="T", "C"="G", "T"="A", "G"="C", "N"="N"))
}

reverseMutation = function(x) {
	.parts = unlist(strsplit(x[["mut"]], split = ""))
	.result = paste(switch(.parts[1], "A"="T", "C"="G", "T"="A", "G"="C", "N":"N"), switch(.parts[2], "A"="T", "C"="G", "T"="A", "G"="C", "N":"N"), sep = "")
	return(.result)
}

data = read_tsv(args[1], col_names = c("pos", "mut", "type", "comptage"))
data$rmut = data$mut
data[data$type == "R", "rmut"] = unlist(apply(data[data$type == "R",], MARGIN = 1, FUN = reverseMutation))

plot = ggplot(data = data[data$pos %in% -50:200,], aes(x=pos, y=comptage, color=type)) + geom_line() + facet_wrap(~rmut) + theme_linedraw() + theme_poster
ggsave("muts_by_side.png", plot, width = 20)

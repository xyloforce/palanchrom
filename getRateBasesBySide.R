library(ggplot2)
library(extrafont)

args = commandArgs(trailingOnly=TRUE)

reverseBase = function(x) {
	return(switch(x[["base"]], "A" = "T", "C" = "G", "T" = "A", "G" = "C", "N" = "N"))
}

source("~/setThemePoster.R")

data = read.delim(args[1], header = F, col.names = c("pos", "base", "type", "comptage"))
data$rbase = data$base
data[data$type == "-", "rbase"] = unlist(apply(data[data$type == "-",], MARGIN = 1, FUN = reverseBase))
head(data)

cumulated = aggregate(data$comptage, by = list(pos = data$pos, type = data$type), FUN = sum)
data = rbind(data, data.frame("pos" = cumulated$pos, "base" = "merged", "type" = cumulated$type, "comptage" = cumulated$x, "rbase" = "merged"))
# data = cumulated

plot = ggplot(data = data[data$pos %in% -50:200,], aes(x = pos, y = comptage, color = rbase)) + geom_line(linewidth = 1) + facet_wrap( ~ type) + theme_poster
ggsave("bases_by_side.png", plot, width = 16, height = 9)

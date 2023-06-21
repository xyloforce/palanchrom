library(stringr)
library(ggplot2)
library(extrafont)

args = commandArgs(trailingOnly=TRUE)
source("~/setThemePoster.R")
# args = c("formattedNCG", "gc_mutations.png")

print("loading df muts")
fileList = list.files(path = args[1], "*_mutations.tsv", full.names = TRUE)

df = data.frame()
for(filename in fileList) {
    df2 = read.delim(filename)
    if(grepl("A|T", filename)) { # want mutations going in the GC direction
        df = rbind(df, cbind(df2$position, "AT", str_match(filename, "([AT])_mutations.tsv")[,2], df2$comptage, rowSums(df2[,grepl("^[AT]", colnames(df2))], na.rm = TRUE)))
    } else { # want mutations going in the AT direction
        df = rbind(df, cbind(df2$position, "GC", str_match(filename, "([GC])_mutations.tsv")[,2], df2$comptage, rowSums(df2[,grepl("^[GC]", colnames(df2))], na.rm = TRUE)))
    }
}
colnames(df) = c("position", "type", "base", "total", "count")

df = aggregate(df[,c("total", "count")], by = list(position = df$position, type = df$type), FUN = function(x) {return(sum(as.numeric(x), na.rm = TRUE))})
df$position = as.numeric(df$position)

df$ratio = df$count / df$total

plot = ggplot(data = df, aes(x = position, y = ratio, color = type)) + geom_line(size = 1) + xlim(-50, 270) + ylim(0.001, 0.01) + theme_poster
ggsave(args[2], width = 20)

library(zoo)
library(ggh4x)
library(stringr)
library(ggplot2)
library(extrafont)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")
# args = c("formattedNCG", "gc_mutations.png", 20)

bin_size = 20
reverse = FALSE
xlimM = -25
xlimP = 250
ylimP = 0.025
ylimM = -0.0025

comptages = FALSE
vertical_bars = FALSE

print("loading df muts")
values = c("ancestral" = args[1], "SNP" = args[2])

final_df = data.frame()
for(context in 1:2) {
	fileList = list.files(path = values[context], "*_mutations.tsv", full.names = TRUE)
	df = data.frame()
	for (filename in fileList) {
	    print(filename)
	    df2 = read.delim(filename)
	    if (grepl("A_mutations|T_mutations", filename)) { # want mutations going in the GC direction
		tmp = cbind(df2$position, "W → S", str_match(filename, "([AT])_mutations.tsv")[, 2], df2$comptage, rowSums(df2[, grepl("^[AT][GC]$", colnames(df2))], na.rm = TRUE))
		tmp = rbind(tmp, cbind(df2$position, "W → W | S → S", str_match(filename, "([AT])_mutations.tsv")[, 2], df2$comptage, df2[, grepl("^[AT][AT]$", colnames(df2))]))
	    } else { # want mutations going in the AT direction
		tmp = cbind(df2$position, "S → W", str_match(filename, "([GC])_mutations.tsv")[, 2], df2$comptage, rowSums(df2[, grepl("^[GC][AT]$", colnames(df2))], na.rm = TRUE))
		tmp = rbind(tmp, cbind(df2$position, "W → W | S → S", str_match(filename, "([GC])_mutations.tsv")[, 2], df2$comptage, df2[, grepl("^[GC][GC]$", colnames(df2))]))
	    }
	    df = rbind(df, tmp)
	}
	df$type = names(values[context])
	final_df = rbind(final_df, df)
}
colnames(final_df) = c("position", "type", "base", "total", "count", "context")
head(final_df)

final_df = aggregate(final_df[, c("total", "count")], by = list(position = final_df$position, type = final_df$type, context = final_df$context), FUN = function(x) {
    return(sum(as.numeric(x), na.rm = TRUE))
})
final_df$position = as.numeric(final_df$position)
final_df = final_df[final_df$position < xlimP + bin_size & final_df$position > xlimM - bin_size, ]
final_df = final_df[order(final_df$position), ]

print("Binning")
matchs = c("total" = "total_bases", "count" = "slidedCount")
for(context in unique(final_df$context)) {
	for(type in unique(final_df$type)) {
		for(kind in 1:2) {
			final_df[final_df$type == type, matchs[kind]] = rollapply(final_df[final_df$type == type, names(matchs[kind])], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")
		}
	}

}

print("Calculating ratios")

final_df$ratio = final_df$slidedCount / final_df$total_bases

selection = reshape(data = final_df[,c("position", "type", "context", "ratio")], direction = "wide",
		idvar = c("position", "type"), timevar = "context",
		v.names = "ratio")
selection[selection$position == 0,]

selection$final_ratio = selection$ratio.ancestral / selection$ratio.SNP

colors = c("W → S" = "#FF713D", "W → W | S → S" = "#20E8C0", "S → W" = "#9C310B")
ggplot(data = selection, aes(x = position, y = final_ratio, color = type)) + geom_line(linewidth = 1) + facet_wrap(~ type) + theme_poster + ylab("Ancestral rate / SNP rate") + scale_color_manual(values = colors) + theme(legend.position = "None")
ggsave("pic.svg", width = 20, height = 9)

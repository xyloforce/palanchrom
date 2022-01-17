#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
#args = c("data/202201141625panTro5_nCG_formatted/", "CG")

theme_bob = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(size = 2),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                strip.text.x = element_text(size = 16),
                strip.text.y = element_text(size = 16),
                axis.text = element_text(size=14),
                #plot.title = element_text(face = "bold"),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 14),
                plot.title = element_text(size = 20, hjust = 0.5))

#### PLOT SOURCE X DEST X CONTEXT ################################################################################################

## load muts according to wildcard

filelist = list.files(path = args[1], "*_mutations.tsv", full.names = TRUE)
df = read_tsv(filelist[1], show_col_types = FALSE)

for(filepath in filelist[2:length(filelist)]) {
	tmp = read_tsv(filepath, show_col_types = FALSE)
	tmp$position = NULL
	df = cbind(df, tmp)
}

## keep only "mean" ones

df = cbind(df$position, df[, grep("mean10", colnames(df))])

## rename cols to make them readable

col_names_new = sapply(str_split(colnames(df[,2:ncol(df)]), "_"), FUN = function(x) x[2])
colnames(df) = c("position", col_names_new)

## reshape to ggplot them

df = reshape(data = df, idvar = "position", varying = c(colnames(df[,2:ncol(df)])), v.name = c("mean10"), times = c(colnames(df[,2:ncol(df)])), direction = "long")
colnames(df) = c("position", "mutation", "mean10")
df = cbind(df, str_split_fixed(df$mutation, pattern ="", n=2))
colnames(df) = c("position", "mutation", "mean10", "ancestral", "reference")

df$group = sapply(df$mutation, FUN = function(x) switch(x, "AT" = 2, "TA" = 2, "AG" = 1, "TC" = 1, "AC" = 3, "TG" = 3, "GC" = 6, "CG" = 6, "GT" = 5, "CA" = 5, "GA" = 4, "CT" = 4))

correct_labels <- c("3" = "AC - TG (transversion)",
                    "2" = "AG - TC (transition)",
                    "1" = "AT - TA (transversion)",
                    "5" = "CA - GT (transversion)",
                    "4" = "CG - GC (transversion)",
                    "6" = "CT - GA (transition)"
                    )

df$color = sapply(df$mutation, FUN = function(x) switch(x, "AT" = "1", "TA" = "2", "AG" = "3", "TC" = "4", "AC" = "1", "TG" = "2", "GC" = "2", "CG" = "1", "GT" = "2", "CA" = "1", "GA" = "4", "CT" = "3"))

# limits = c(0.0001, 0.01), breaks = seq(0.0001, 0.01, by = 0.0015),
# how to add breaks AND scales free ?
plot1 = ggplot(data = df[df$position %in% -50:500,], aes(y = mean10, x = position, color = color)) +
	facet_wrap(~group, labeller = as_labeller(correct_labels), scales = "free") + geom_line() +
	ylab("% de mutation lissés sur 10 pb") +
	scale_color_discrete(type = c("#7F0794", "#E033FF", "#509400", "#84E016")) +
	ggtitle("Les taux de transversion sont supérieurs aux taux de transition") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	geom_vline(xintercept = c(0, 117, 270), color = "grey") +
	theme_linedraw() +
	theme(strip.placement = "outside", legend.position = "none") +
	theme_bob

## load data to add second plot of global muts

filepath = list.files(path = args[1], "*_total.tsv", full.names = TRUE)

total = read_tsv(filepath[1], show_col_types = FALSE)

plot2 = ggplot(data = total[total$position %in% -50:500,], aes(x=position, y = mean10)) + geom_line() +
	ylab("% de mutation lissés sur 10 pb") +
	ggtitle("Taux de mutation global") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	geom_vline(xintercept = c(0, 117, 270), color = "grey") +
	theme_linedraw() +
	theme_bob

figure = plot_grid(plot1, plot2, labels = c("A", "B"), rel_widths = c(0.6, 0.3))
figure

folder = args[2]

if(file.exists(folder)) {
	print("Folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

ggsave("bases_by_group.png", figure, width = 24, height = 9)

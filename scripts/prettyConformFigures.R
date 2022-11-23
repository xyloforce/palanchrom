#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)
library(extrafont)

args = commandArgs(trailingOnly=TRUE)
#args = c("formatted", "figure")

source("~/setThemePoster.R")
# theme_bob = theme(panel.border = element_rect(size = 2),
#                 axis.title.x = element_text(size = 16),
#                 axis.title.y = element_text(size = 16),
#                 strip.text.x = element_text(size = 16),
#                 strip.text.y = element_text(size = 16),
#                 axis.text = element_text(size=14),
#                 #plot.title = element_text(face = "bold"),
#                 legend.title = element_text(size = 16),
#                 legend.text = element_text(size = 14),
#                 plot.title = element_text(size = 20, hjust = 0.5))

#### PLOT SOURCE X DEST X CONTEXT ################################################################################################

## load muts according to wildcard

print("loading df muts")
filelist = list.files(path = args[1], "*_mutations.tsv", full.names = TRUE)
df_source = read_tsv(filelist[1], show_col_types = FALSE)

if(length(args) > 2) {
	xlim1 = as.numeric(args[3])
	xlim2 = as.numeric(args[4])
} else {
	xlim1 = -50
	xlim2 = 350
}

for(filepath in filelist[2:length(filelist)]) {
	tmp = read_tsv(filepath, show_col_types = FALSE)
	tmp$position = NULL
	df_source = cbind(df_source, tmp)
}

## keep only "mean" ones

print("selecting columns")
df = cbind(df_source$position, df_source[, grep("mean10", colnames(df_source))])
#df2 = cbind(df_source$position, df_source[, grep("error10", colnames(df_source))])
df2 = cbind(df_source$position, df_source[, grep("error10", colnames(df_source))])
## rename cols to make them readable

print("renaming cols")
col_names_new = sapply(str_split(colnames(df[,2:ncol(df)]), "_"), FUN = function(x) x[2])
colnames(df) = c("position", col_names_new)
col_names_new = sapply(str_split(colnames(df2[,2:ncol(df2)]), "_"), FUN = function(x) x[2])
colnames(df2) = c("position", col_names_new)
df2 = df2[df2$position %% 10 == 0,]

## reshape to ggplot them

print("reshaping")
df = reshape(data = df, idvar = "position", varying = c(colnames(df[,2:ncol(df)])), v.name = c("mean10"), times = c(colnames(df[,2:ncol(df)])), direction = "long")
colnames(df) = c("position", "mutation", "mean10")
df = cbind(df, str_split_fixed(df$mutation, pattern ="", n=2))
colnames(df) = c("position", "mutation", "mean10", "ancestral", "reference")
df2 = reshape(data = df2, idvar = "position", varying = c(colnames(df2[,2:ncol(df2)])), v.name = c("error10"), times = c(colnames(df2[,2:ncol(df2)])), direction = "long")
colnames(df2) = c("position", "mutation", "error10")
df$error10 = df2[match(df$position, df2$position), "error10"]
df$ymax = df$mean10 + df$error10
df$ymin = df$mean10 - df$error10

df$group = sapply(df$mutation, FUN = function(x) switch(x, "AT" = 2, "TA" = 2, "AG" = 1, "TC" = 1, "AC" = 3, "TG" = 3, "GC" = 6, "CG" = 6, "GT" = 5, "CA" = 5, "GA" = 4, "CT" = 4))

correct_labels <- c("3" = "A→C - T→G (transversion)",
                    "1" = "A→G - T→C (transition)",
                    "2" = "A→T - T→A (transversion)",
                    "5" = "C→A - G→T (transversion)",
                    "6" = "C→G - G→C (transversion)",
                    "4" = "C→T - G→A (transition)"
                    )

df$color = sapply(df$mutation, FUN = function(x) switch(x, "AT" = "1", "TA" = "2", "AG" = "3", "TC" = "4", "AC" = "1", "TG" = "2", "GC" = "2", "CG" = "1", "GT" = "2", "CA" = "1", "GA" = "4", "CT" = "3"))

print("plotting")
df = df[df$position > xlim1 & df$position < xlim2,]
# limits = c(0.0001, 0.01), breaks = seq(0.0001, 0.01, by = 0.0015),
# how to add breaks AND scales free ?
plot1 = ggplot(data = df, aes(y = mean10, x = position, color = color)) +
	facet_wrap(~group, labeller = as_labeller(correct_labels), scales = "free") + geom_line(size = 1) +
	geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
	ylab("% de mutation lissés sur 10 pb") +
	scale_color_discrete(type = c("#9C310B", "#FF713D", "#009C7D", "#20E8C0")) +
	ggtitle("Taux de mutation complémentaires (le premier est le plus sombre)") +
	scale_x_continuous(limits = c(xlim1, xlim2), sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	geom_vline(xintercept = c(0, 133, 266), color = "grey") +
	theme_poster +
	theme(strip.placement = "outside", legend.position = "none")

## load data to add second plot of global muts

filepath = list.files(path = args[1], "total.tsv", full.names = TRUE)

total = read_tsv(filepath[1], show_col_types = FALSE)

total = total[total$position > xlim1 & total$position < xlim2,]

plot2 = ggplot(data = total, aes(x=position, y = mean10)) + geom_line(size = 1) +
	ylab("% de mutation lissés sur 10 pb") +
	ggtitle("Taux de mutation global") +
	scale_x_continuous(limits = c(xlim1, xlim2), sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	geom_vline(xintercept = c(0, 133, 266), color = "grey") +
	theme_poster

figure = plot_grid(plot1, plot2, labels = c("A", "B"), rel_widths = c(0.6, 0.3))
# figure

folder = args[2]

if(file.exists(folder)) {
	print("Folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

ggsave("bases_by_group.png", figure, width = 24, height = 9)

#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)

#args = commandArgs(trailingOnly=TRUE)
args = c("hg38.rda", "panTro5.rda")

filelist = list.files(path = args[1], "*_total.tsv", full.names = TRUE)
df = read_tsv(filelist[1])
species = str_split(args[1], "/")
species = str_split(species[length(species)], ".")[1]
df$species = species

filelist = list.files(path = args[2], "*_total.tsv", full.names = TRUE)
tmp = read_tsv(filelist[1])
species = str_split(args[2], "/")
species = str_split(species[length(species)], ".")[1]
tmp$species = species

df = rbind(df, tmp)

plot = ggplot(data = total[total$position %in% -50:500,], aes(x=position, y = mean10, color = species)) + geom_line() +
	ylab("% de mutation lissés sur 10 pb") +
	ggtitle("Taux de mutation global") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	theme_linedraw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = c(0, 117, 270), color = "grey")

folder = args[3]

if(file.exists(folder)) {
	print("Folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

ggsave("mutations_by_species.png", figure, width = 24, height = 9)

#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)

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
                legend.text = element_text(size = 14))
                #, legend.position = c(0.8, 0.3)

args = commandArgs(trailingOnly=TRUE)
#args = c("data/202201171632_panTro5_CG_formatted", "data/202201171632_hg38_CG_formatted", "data/202201171632_panTro5_hg38_CG_figs")

filelist = list.files(path = args[1], "*_total.tsv", full.names = TRUE)
print(filelist)
df = read_tsv(filelist[1], show_col_types = FALSE)
print(head(df))
species = unlist(str_split(args[1], "/"))
species = unlist(str_split(species[length(species)], "_"))[2]
df$species = species

filelist = list.files(path = args[2], "*_total.tsv", full.names = TRUE)
tmp = read_tsv(filelist[1], show_col_types = FALSE)
species = unlist(str_split(args[2], "/"))
species = unlist(str_split(species[length(species)], "_"))[2]
tmp$species = species

df = rbind(df, tmp)
print(head(df))
plot1 = ggplot(data = df[df$position %in% -50:500,], aes(x=position, y = mean10, color = species)) + geom_line() +
	ylab("% de mutation lissés sur 10 pb") +
	ggtitle("Taux de mutation global") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	theme_linedraw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = c(0, 117, 270), color = "grey") + theme_bob
plot1

folder = args[3]

if(file.exists(folder)) {
	print("Folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

ggsave("mutations_by_species.png", plot1, width = 24, height = 9)

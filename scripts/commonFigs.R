#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
#args = c("hg38.rda", "panTro5.rda")

load(args[1])

muts_hg = muts
bases_hg = bases

muts_hg$species = "hg38"
bases_hg$species = "hg38"

load(args[2])

muts_pan = muts
bases_pan = bases

muts_pan$species = "panTro5"
bases_pan$species = "panTro5"

muts = rbind(muts_hg, muts_pan)
muts = aggregate(muts$mean10, by = list(muts$position, muts$species), FUN = sum)
colnames(muts) = c("position", "species", "mean10")

ggplot(data = muts, aes(x = position, y = mean10, color = species)) +
	geom_line() +
	ylab("% de mutation lissés sur 10 pb") +
	ggtitle("Taux de mutation global par espèce") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	theme_linedraw() +
	theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) +
	geom_vline(xintercept = c(0, 117, 270), color = "grey")

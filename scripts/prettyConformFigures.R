#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
#args = c("data/hg38_CPG_bases.tsv", "data/hg38_nCPG_bases.tsv", "data/hg38_CPG_muts.tsv", "data/hg38_nCPG_muts.tsv", "results_temp", "hg38.rda")

mean10pb <- function(x, n = 10){if(length(x) > 0) {return(filter(x, rep(1 / n, n), sides = 2, circular = FALSE))} else {return(NA)}}
print("loading files")
basesCPG = read_tsv(args[1], col_names = c("position", "base", "type", "comptage"), na = character())
basesNCPG = read_tsv(args[2], col_names = c("position", "base", "type", "comptage"), na = character())
mutsCPG = read_tsv(args[3], col_names = c("position", "mutation", "type", "comptage"), na = character())
mutsNCPG = read_tsv(args[4], col_names = c("position", "mutation", "type", "comptage"), na = character())

folder = args[5]

if(file.exists(folder)) {
	print("Folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

mutsCPG$source = "CPG"
mutsNCPG$source = "NCPG"
basesCPG$source = "CPG"
basesNCPG$source = "NCPG"

### CHECK THAT L = R FROM THE MUTATION RATE POINT OF VIEW
print("plotting with type")
muts = rbind(mutsCPG, mutsNCPG)
muts = cbind(muts, str_split_fixed(muts$mutation, "", n=2))
colnames(muts) = c("position", "mutation", "type", "comptage", "source", "bsource", "bdestination")
bases = rbind(basesCPG, basesNCPG)

muts = muts[muts$position %in% -90:650,]
muts$indexes = 0
for(base in unique(muts$bsource)) {
	for(type in unique(muts$type)) {
		for(source in unique(muts$source)) {
			for(position in unique(muts$position)) {
				muts[muts$position == position & muts$bsource == base & muts$source == source & muts$type == type,"indexes"] = which(bases$position == position & bases$source == source & bases$base == base & bases$type == type)
			}
		}
	}
}

muts = as.data.frame(muts)
muts$relative = 0
muts$relative = unlist(muts$comptage / bases[muts$indexes,"comptage"] * 100)

print("averaging results to soften the curve")
muts$mean10 = 0
for(base in unique(muts$bsource)) {
	for(type in unique(muts$type)) {
		for(source in unique(muts$source)) {
			muts[muts$bsource == base & muts$type == type & muts$source == source, "mean10"] = mean10pb(muts[muts$bsource == base & muts$type == type & muts$source == source, "relative"])
		}
	}
}
muts = muts[muts$position %in% -80:550,]

ggplot(data = muts[muts$source == "NCPG",], aes(x = position, y = mean10, color = type)) +
	geom_line() + facet_wrap(bsource ~ bdestination) +
	ggtitle("Mutations (nCPG) en fonction du côté\nde l'intervalle considéré") +
	ylab("% de mutation lissés sur 10 pb") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	theme_linedraw() + theme(strip.placement = "outside", panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = c(0, 117, 270), color = "grey")

ggsave("mutations rates by side.png", width = 10, height = 10)
### CHECK THAT GC MUTS = CG MUTS GLOBALLY
print("merging types with reverse complement")
reverseBase = function(x) {
	return(switch(x[["base"]], "A"="T", "C"="G", "T"="A", "G"="C", "N"="N"))
}
reverseMutation = function(x) {
	.parts = unlist(strsplit(x[["mutation"]], split = ""))
	.result = paste(switch(.parts[1], "A"="T", "C"="G", "T"="A", "G"="C", "N":"N"), switch(.parts[2], "A"="T", "C"="G", "T"="A", "G"="C", "N":"N"), sep = "")
	return(.result)
}

muts = rbind(mutsCPG, mutsNCPG)
bases = rbind(basesCPG, basesNCPG)
muts[muts$type == "R","mutation"] = unlist(apply(muts[muts$type == "R",], MARGIN = 1, FUN = reverseMutation))
muts = aggregate(muts$comptage, by = list(muts$position, muts$mutation, muts$source), FUN = sum)
colnames(muts) = c("position", "mutation", "source", "comptage")
bases[bases$type == "R","base"] = unlist(apply(bases[bases$type == "R",], MARGIN = 1, FUN = reverseBase))
bases = aggregate(bases$comptage, by = list(bases$position, bases$base, bases$source), FUN = sum)
colnames(bases) = c("position", "base", "source", "comptage")

muts = cbind(muts, str_split_fixed(muts$mutation, "", n=2))
colnames(muts) = c("position", "mutation", "source", "comptage", "bsource", "bdestination")

muts = muts[muts$position %in% -90:650,]
muts$indexes = 0
for(base in unique(muts$bsource)) {
	for(source in unique(muts$source)) {
		for(position in unique(muts$position)) {
			muts[muts$position == position & muts$bsource == base & muts$source == source,"indexes"] = which(bases$position == position & bases$source == source & bases$base == base)
		}
	}
}

muts$relative = muts$comptage / bases[muts$indexes,"comptage"] * 100

for(base in unique(muts$bsource)) {
	for(source in unique(muts$source)) {
		muts[muts$bsource == base & muts$source == source, "mean10"] = mean10pb(muts[muts$bsource == base & muts$source == source, "relative"])
	}
}
muts = muts[muts$position %in% -80:550,]

muts$group = sapply(muts$mutation, FUN = function(x) switch(x, "AT" = 1, "TA" = 1, "AG" = 2, "TC" = 2, "AC" = 3, "TG" = 3, "GC" = 4, "CG" = 4, "GT" = 5, "CA" = 5, "GA" = 6, "CT" = 6))

correct_labels <- c("3" = "AC - TG (transversion)",
                    "2" = "AG - TC (transition)",
                    "1" = "AT - TA (transversion)",
                    "5" = "CA - GT (transversion)",
                    "4" = "CG - GC (transversion)",
                    "6" = "CT - GA (transition)"
                    )

plot1 = ggplot(data = muts[muts$source == "NCPG",], aes(y = mean10, x = position, color = bdestination)) +
	facet_wrap(~group, labeller = as_labeller(correct_labels)) + geom_line() +
	geom_point(data = muts[muts$source == "NCPG" & muts$position %% 30 == 0,], aes(y = mean10, x = position, shape = bsource)) +
	ylab("% de mutation lissés sur 10 pb") +
	scale_color_discrete(name = "Base du génome\nde référence") +
	scale_shape_discrete(name = "Base du génome\nancestral") +
	ggtitle("Les taux de transversion sont supérieurs aux taux de transition") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	theme_linedraw() + theme(strip.placement = "outside", panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = c(0, 117, 270), color = "grey")

tmp = aggregate(muts$mean10, by = list(muts$position, muts$source), FUN = sum)
colnames(tmp) = c("position", "source", "mean10")

plot2 = ggplot(data = tmp[tmp$source == "NCPG",], aes(y = mean10, x = position)) +
	geom_line() +
	ylab("% de mutation lissés sur 10 pb") +
	ggtitle("Taux de mutation global") +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	theme_linedraw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = c(0, 117, 270), color = "grey")

#figure <- ggdraw() + draw_plot(plot2, 0, 0, 1, 1) + draw_plot(plot1, 1, 0, 2, 6) + draw_plot_label(c("A", "B"), c(0, 0), c(0, 1), size = 15)
figure = plot_grid(plot1, plot2, labels = c("A", "B"), rel_widths = c(1, 0.3))
figure
ggsave("mutations rates by base and global.png", plot = figure, width = 15, height = 5)

save(muts, bases, file = args[6])

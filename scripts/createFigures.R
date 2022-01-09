#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

mean10pb <- function(x, n = 10){if(length(x) > 0) {return(filter(x, rep(1 / n, n), sides = 2, circular = FALSE))} else {return(NA)}}
print("loading files")
basesCPG = read_tsv(args[1], col_names = c("position", "base", "type", "comptage"), na = character())
basesNCPG = read_tsv(args[2], col_names = c("position", "base", "type", "comptage"), na = character())
mutsCPG = read_tsv(args[3], col_names = c("position", "mutation", "type", "comptage"), na = character())
mutsNCPG = read_tsv(args[4], col_names = c("position", "mutation", "type", "comptage"), na = character())

folder = args[5]

if(file.exists(folder)) {
	stop("Error : folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

mutsCPG$source = "CPG"
mutsNCPG$source = "NCPG"
basesCPG$source = "CPG"
basesNCPG$source = "NCPG"

############# TEST WITH TYPE : GET RATE : BASE / COUNT FOR BASE #############
print("plotting with type")
muts = rbind(mutsCPG, mutsNCPG)
muts = cbind(muts, str_split_fixed(muts$mutation, "", n=2))
colnames(muts) = c("position", "mutation", "type", "comptage", "source", "bsource", "bdestination")
head(muts)
bases = rbind(basesCPG, basesNCPG)

muts = muts[muts$position %in% -90:650,]
muts$indexes = 0
print(unique(muts$bsource))
print(unique(muts$mutation))
print(muts[is.na(muts$mutation),])
print(bases[bases$position == -90 & bases$type == "L" & bases$source == "NCPG",])
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
			print(muts[muts$bsource == base & muts$type == type & muts$source == source, "relative"])
			muts[muts$bsource == base & muts$type == type & muts$source == source, "mean10"] = mean10pb(muts[muts$bsource == base & muts$type == type & muts$source == source, "relative"])
		}
	}
}
muts = muts[muts$position %in% -80:550,]

ggplot(data = muts[muts$source == "NCPG",], aes(x = position, y = mean10, color = type)) + geom_line() + facet_wrap(bsource ~ bdestination) + theme_minimal() + theme(panel.grid = element_line(linetype = "solid", size = 0.5, color = "#868686"))
ggsave("relative_muts_by_source_dest_type.png", width = 10, height = 10, bg = "white")

############# GET INFO ABOUT COUNTS #############
print("plotting counts")
tmp = rbind(basesCPG, basesNCPG)
tmp = aggregate(tmp$comptage, by=list(tmp$position, tmp$source), FUN = sum)
colnames(tmp) = c("position", "source", "comptage")
ggplot(data = tmp[tmp$position %in% -50:500,], aes(x=position, y = comptage, color = source)) + geom_point() + scale_y_log10()  + theme_minimal() + theme(panel.grid = element_line(linetype = "solid", size = 0.5, color = "#868686"))
ggsave("global_counts_by_pos.png", width = 10, height = 5, bg = "white")

############# REVERT BASE TO DELETE TYPE INFO #############
print("merging types with reverse complement")
reverseBase = function(x) {
	return(switch(x[["base"]], "A"="T", "C"="G", "T"="A", "G"="C", "N"="N"))
}
basesCPG[basesCPG$type == "R","base"] = unlist(apply(basesCPG[basesCPG$type == "R",], MARGIN = 1, FUN = reverseBase))
basesCPG = aggregate(basesCPG$comptage, by = list(basesCPG$position, basesCPG$base, basesCPG$source), FUN = sum)
colnames(basesCPG) = c("position", "base", "source", "comptage")

basesNCPG[basesNCPG$type == "R","base"] = unlist(apply(basesNCPG[basesNCPG$type == "R",], MARGIN = 1, FUN = reverseBase))
basesNCPG = aggregate(basesNCPG$comptage, by = list(basesNCPG$position, basesNCPG$base, basesNCPG$source), FUN = sum)
colnames(basesNCPG) = c("position", "base", "source", "comptage")

############# REVERT MUTS TO DELETE TYPE INFO #############
reverseMutation = function(x) {
	.parts = unlist(strsplit(x[["mutation"]], split = ""))
	.result = paste(switch(.parts[1], "A"="T", "C"="G", "T"="A", "G"="C", "N":"N"), switch(.parts[2], "A"="T", "C"="G", "T"="A", "G"="C", "N":"N"), sep = "")
	return(.result)
}

mutsNCPG[mutsNCPG$type == "R","mutation"] = unlist(apply(mutsNCPG[mutsNCPG$type == "R",], MARGIN = 1, FUN = reverseMutation))
mutsNCPG = aggregate(mutsNCPG$comptage, by = list(mutsNCPG$position, mutsNCPG$mutation, mutsNCPG$source), FUN = sum)
colnames(mutsNCPG) = c("position", "mutation", "source", "comptage")

mutsCPG[mutsCPG$type == "R","mutation"] = unlist(apply(mutsCPG[mutsCPG$type == "R",], MARGIN = 1, FUN = reverseMutation))
mutsCPG = aggregate(mutsCPG$comptage, by = list(mutsCPG$position, mutsCPG$mutation, mutsCPG$source), FUN = sum)
colnames(mutsCPG) = c("position", "mutation", "source", "comptage")

############# SEPARATE BASES IN MUT #############
mutsCPG = cbind(mutsCPG, str_split_fixed(mutsCPG$mutation, "", n=2))
colnames(mutsCPG) = c("position", "mutation", "source", "comptage", "bsource", "bdestination")
mutsNCPG = cbind(mutsNCPG, str_split_fixed(mutsNCPG$mutation, "", n=2))
colnames(mutsNCPG) = c("position", "mutation", "source", "comptage", "bsource", "bdestination")

############# GET RATE : BASE / COUNT FOR BASE #############
muts = rbind(mutsCPG, mutsNCPG)
bases = rbind(basesCPG, basesNCPG)

muts = muts[muts$position %in% -90:650,]
muts$indexes = 0
for(base in unique(muts$bsource)) {
	for(source in unique(muts$source)) {
		for(position in unique(muts$position)) {
			muts[muts$position == position & muts$bsource == base & muts$source == source,"indexes"] = which(bases$position == position & bases$source == source & bases$base == base)
		}
	}
}

muts$relative = 0
muts$relative = muts$comptage / bases[muts$indexes,"comptage"] * 100

############# SHOW MUTS BY SOURCE AND DEST #############
print("plotting muts")
for(base in unique(muts$bsource)) {
	for(source in unique(muts$source)) {
		print(muts[muts$bsource == base & muts$source == source, "relative"])
		muts[muts$bsource == base & muts$source == source, "mean10"] = mean10pb(muts[muts$bsource == base & muts$source == source, "relative"])
	}
}
muts = muts[muts$position %in% -80:550,]

ggplot(data = muts, aes(x = position, y = mean10, color = bdestination)) + geom_line() + facet_wrap(source~bsource, scales = "free")  + theme_minimal() + theme(panel.grid = element_line(linetype = "solid", size = 0.5, color = "#868686"))
ggsave("relative_counts_by_pos_and_type.png", width = 30, height = 15, bg = "white")

ggplot(data = muts[muts$source == "CPG",], aes(x = position, y = mean10, color = bdestination)) + geom_line() + facet_wrap(~bsource) + theme_minimal() + theme(panel.grid = element_line(linetype = "solid", size = 0.5, color = "#868686"))
ggsave("relative_counts_by_pos_and_type_CPG.png", width = 30, height = 15, bg = "white")

ggplot(data = muts[muts$source == "NCPG",], aes(x = position, y = mean10, color = bdestination)) + geom_line() + facet_wrap(~bsource)  + theme_minimal() + theme(panel.grid = element_line(linetype = "solid", size = 0.5, color = "#868686"))
ggsave("relative_counts_by_pos_and_type_NCPG.png", width = 10, height = 10, bg = "white")


############# SHOW MUTS BY SOURCE #############
print("plotting muts for bsource only")
tmp = aggregate(muts$relative, by=list(muts$position, muts$source, muts$bsource), FUN = sum)
colnames(tmp) = c("position", "source", "bsource", "relative")
ggplot(data= tmp, aes(x=position, y=relative, color = bsource)) + geom_point() + facet_wrap(~source, scales = "free") + theme_minimal()
ggsave("relative_counts_by_pos_and_base.png", width = 30, height = 15, bg = "white")

############# SHOW MUTS BY POS #############
########### FIRST RECALCUlATE % ############
print("plotting mutation rates")
muts = rbind(mutsCPG, mutsNCPG)
bases = rbind(basesCPG, basesNCPG)

muts = muts[muts$source == "NCPG",]
muts = aggregate(muts$comptage, by=list(muts$position), FUN = sum)
colnames(muts) = c("position", "comptage")
bases = bases[bases$source == "NCPG",]
bases = aggregate(bases$comptage, by=list(bases$position), FUN = sum)
colnames(bases) = c("position", "comptage")

muts = muts[muts$position %in% -50:500,]
muts$indexes = 0
muts$indexes = match(muts$position, bases$position)
muts$relative = 0
muts$relative = muts$comptage / bases[muts$indexes,"comptage"] * 100
muts$mean10 = mean10pb(muts$relative)

############### THEN PLOT ###############
ggplot(data= muts, aes(x=position, y=mean10)) + geom_line() + theme_minimal()
ggsave("relative_counts_by_pos.png", width = 30, height = 30, bg = "white")

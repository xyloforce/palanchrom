#!/usr/bin/env Rscript
library(readr)
library(ggplot2)
library(stringr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
#args = c("data/hg38_nCPG_bases.tsv", "data/hg38_nCPG_muts.tsv")

mean10pb <- function(x, n = 10){if(length(x) > 0) {return(filter(x, rep(1 / n, n), sides = 2, circular = FALSE))} else {return(NA)}}

print("loading files")
bases = read_tsv(args[1], col_names = c("position", "base", "type", "comptage"), na = c(""), show_col_types = FALSE)
muts = read_tsv(args[2], col_names = c("position", "mutation", "type", "comptage"), na = c(""), show_col_types = FALSE)

folder = args[3]
if(file.exists(folder)) {
	print("Folder exists")
} else {
	dir.create(folder)
}

setwd(folder)

bases = bases[bases$position <= 5005,]
muts = muts[muts$position <= 5005,]

#### FIRST : REVERSE TO DELETE TYPE INFO ########################################################################

reverseBase = function(x) {
	return(switch(x[["base"]], "A"="T", "C"="G", "T"="A", "G"="C", "N"="N"))
}
reverseMutation = function(x) {
	.parts = unlist(strsplit(x[["mutation"]], split = ""))
	.result = paste(switch(.parts[1], "A"="T", "C"="G", "T"="A", "G"="C", "N"="N"), switch(.parts[2], "A"="T", "C"="G", "T"="A", "G"="C", "N"="N"), sep = "")
	return(.result)
}

print("get all bases on one side")
muts[muts$type == "R","mutation"] = unlist(apply(muts[muts$type == "R",], MARGIN = 1, FUN = reverseMutation))
muts = aggregate(muts$comptage, by = list(muts$position, muts$mutation), FUN = sum)
colnames(muts) = c("position", "mutation", "comptage")
bases[bases$type == "R","base"] = unlist(apply(bases[bases$type == "R",], MARGIN = 1, FUN = reverseBase))
bases = aggregate(bases$comptage, by = list(bases$position, bases$base), FUN = sum)
colnames(bases) = c("position", "base", "comptage")

#### MERGE TO CREATE TOTAL DF ########################################################################

print("creating total df")
total = aggregate(bases$comptage, by = list(bases$position), FUN = sum)
colnames(total) = c("position", "comptage")
tmp = aggregate(muts$comptage, by = list(muts$position), FUN = sum)
colnames(tmp) = c("position", "comptage")
total$mutations = tmp[match(total$position, tmp$position),"comptage"]
head(total)
if (unique(-225:5005 %in% total$position)) {
	print("ok")
} else {
	tmp = data.frame(position = -225:5005, comptage=0, mutations=0)
	tmp[match(total$position, tmp$position),] = total
	total = tmp
	total[is.na(total)] = 0
} # prepare code for CPGs

total$relative = total$mutations/total$comptage * 100
total$mean10 = mean10pb(total$relative)
write_tsv(total[total$position %in% -220:5000,], paste("total.tsv", sep = ""))

#### CREATE ONE DF BY TYPE OF MUT ########################################################################
print("create one df by type of mut")
muts = cbind(muts, str_split_fixed(muts$mutation, "", n=2))
colnames(muts) = c("position", "mutation", "comptage", "ancestral", "reference")
for(base in unique(muts$ancestral)) {
	print(base)
	currentBDF = data.frame(position = unique(total$position))
	currentBDF[match(bases[bases$base == base, "position"], currentBDF$position), "comptage"] = bases[bases$base == base, "comptage"]
	for(bdest in unique(muts[muts$ancestral == base, "reference"])) {
		mut = paste(base, bdest, sep = "")
		currentBDF[match(muts[muts$mutation == mut, "position"], currentBDF$position),mut] = muts[muts$mutation == mut, "comptage"]
		currentBDF[is.na(currentBDF[,mut]),mut] = 0
		relative = paste("relative_", mut, sep= "")
		currentBDF[,relative] = currentBDF[,mut] / currentBDF$comptage * 100
		mean10m = paste("mean10_", mut, sep = "")
		currentBDF[,mean10m] = mean10pb(currentBDF[,relative])
	}
	write_tsv(currentBDF[currentBDF$position %in% -220:5000,], paste(base, "_mutations.tsv", sep = ""))
}

#### CREATE ONE DF BY TYPE OF COMPL MUTS ########################################################################
print("create one df by type of complementary mut")
muts = muts[!(grepl("N", muts$mutation)),]
muts$group = sapply(muts$mutation, FUN = function(x) switch(x, "AT" = 1, "TA" = 1, "AG" = 2, "TC" = 2, "AC" = 3, "TG" = 3, "GC" = 4, "CG" = 4, "GT" = 5, "CA" = 5, "GA" = 6, "CT" = 6))

for(group in unique(muts$group)) {
	currentBDF = data.frame(position = unique(total$position))
	print(group)
	for(base in unique(muts[muts$group == group,"ancestral"])) {
		c_base = paste("comptage_", base, sep = "")
		currentBDF[match(bases[bases$base == base, "position"], currentBDF$position), c_base] = bases[bases$base == base, "comptage"]
		for(bdest in unique(muts[muts$ancestral == base & muts$group == group, "reference"])) {
			mut = paste(base, bdest, sep = "")
			currentBDF[match(muts[muts$mutation == mut, "position"], currentBDF$position),mut] = muts[muts$mutation == mut, "comptage"]
			currentBDF[is.na(currentBDF[,mut]),mut] = 0
		}
	}
	print(head(currentBDF))
	currentBDF$relative_both = (currentBDF[,3] + currentBDF[,5]) / (currentBDF[,2] + currentBDF[,4]) * 100
	currentBDF$mean10 = mean10pb(currentBDF$relative_both)
	write_tsv(currentBDF[currentBDF$position %in% -220:5000,], paste("group", group, "-", mut, "_and_reverse.tsv", sep = ""))
}

library(zoo)
library(readr)
library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
#args = c("basesnCG.tsv", "mutsnCG.tsv", "formattednCG/")

#mean10pb <- function(x, n = 10){if(length(x) > 0) {return(filter(x, rep(1 / n, n), sides = 2, circular = FALSE))} else {return(NA)}}
mean10pb = function(x, n = 10) {
    rollapply(x, n, mean, na.rm = TRUE, fill = NA)
}

print("loading files")
bases = read_tsv(args[1], col_names = c("position", "base", "type", "comptage"), na = c(""), col_types = "icci")
muts = read_tsv(args[2], col_names = c("position", "mutation", "type", "comptage"), na = c(""), show_col_types = FALSE)

folder = args[3]
if (file.exists(folder)) {
    print("Folder exists")
} else {
    dir.create(folder)
}

setwd(folder)

bases = bases[bases$position <= 5005, ]
muts = muts[muts$position <= 5005, ]

#### FIRST : REVERSE TO DELETE TYPE INFO ########################################################################

reverse_base = function(col_base) {
    result = sapply(col_base, function(x) {
        switch(x, "A" = "T", "C" = "G", "T" = "A", "G" = "C", "N" = "N")
    })
    return(result)
}
reverse_mutation = function(col_mut) {
    .parts = unlist(str_split(col_mut, pattern = "", simplify = TRUE))
    result = paste(reverse_base(.parts[, 1]), reverse_base(.parts[, 2]), sep = "")
    return(result)
}

if (nrow(muts[muts$type == "-", ]) > 0) {
    print("get all bases on one side")
    muts[muts$type == "-", "mutation"] = sapply(muts[muts$type == "-", "mutation"], FUN = reverse_mutation)
    bases[bases$type == "-", "base"] = sapply(bases[bases$type == "-", "base"], FUN = reverse_base)
    bases = aggregate(bases$comptage, by = list(bases$position, bases$base), FUN = sum, na.rm = TRUE)
    colnames(bases) = c("position", "base", "comptage")
} else {
    print("only one side provided")
}
muts = aggregate(muts$comptage, by = list(muts$position, muts$mutation), FUN = sum, na.rm = TRUE)
colnames(muts) = c("position", "mutation", "comptage")
muts = muts[muts$position > -225 & muts$position < 5005, ]
bases = aggregate(bases$comptage, by = list(bases$position, bases$base), FUN = sum, na.rm = TRUE)
colnames(bases) = c("position", "base", "comptage")
bases = bases[bases$position > -225 & bases$position < 5005, ] 
#### MERGE TO CREATE TOTAL DF ########################################################################

print("creating total df")
total = aggregate(bases$comptage, by = list(bases$position), FUN = sum, na.rm = TRUE)
colnames(total) = c("position", "comptage")
tmp = aggregate(muts$comptage, by = list(muts$position), FUN = sum, na.rm = TRUE)
colnames(tmp) = c("position", "comptage")
total$mutations = 0
total$mutations = tmp[match(total$position, tmp$position), "comptage"]
# total[is.na(total)] = 0
head(total)
if (unique(-225:5005 %in% total$position)) {
    print("ok")
} else {
    tmp = data.frame(position = -225:5005, comptage = 0, mutations = 0)
    tmp[match(total$position, tmp$position), ] = total
#     tmp[is.na(tmp)] = 0
    total = tmp
} # prepare code for CPGs

total$relative = total$mutations / total$comptage
total$error_bar = (2 * sqrt(total$relative * total$comptage * (1 - total$relative))) / total$comptage
total$error10 = mean10pb(total$error_bar)
total$mean10 = mean10pb(total$relative)
write_tsv(total[total$position %in% -220:5000, ], paste("total.tsv", sep = ""))

#### CREATE ONE DF BY TYPE OF MUT ########################################################################
print("create one df by type of mut")
muts = cbind(muts, str_split_fixed(muts$mutation, "", n = 2))
colnames(muts) = c("position", "mutation", "comptage", "ancestral", "reference")
for (base in unique(muts$ancestral)) {
    print(base)
    current_base_df = data.frame(position = unique(total$position), comptage = 0)
    current_base_df[match(bases[bases$base == base, "position"], current_base_df$position), "comptage"] = bases[bases$base == base, "comptage"]
    for (bdest in unique(muts[muts$ancestral == base, "reference"])) {
        mut = paste(base, bdest, sep = "")
        current_base_df[, mut] = 0
        current_base_df[match(muts[muts$mutation == mut, "position"], current_base_df$position), mut] = muts[muts$mutation == mut, "comptage"]
        relative = paste("relative_", mut, sep = "")
        current_base_df[, relative] = current_base_df[, mut] / current_base_df$comptage
        error_bar = paste("error_", mut, sep = "")
        current_base_df[, error_bar] = (2 * sqrt(current_base_df[, relative] * current_base_df$comptage * (1 - current_base_df[, relative]))) / current_base_df$comptage # sqrt(np*(1-p))
        mean10m = paste("mean10_", mut, sep = "")
        current_base_df[, mean10m] = mean10pb(current_base_df[, relative])
        error10m = paste("error10_", mut, sep = "")
        current_base_df[, error10m] = mean10pb(current_base_df[, error_bar])
    }
#     currentBDF[is.na(currentBDF)] = 0
    write_tsv(current_base_df[current_base_df$position %in% -220:5000, ], paste(base, "_mutations.tsv", sep = ""))
}

#### CREATE ONE DF BY TYPE OF COMPL MUTS ########################################################################
print("create one df by type of complementary mut")
muts = muts[!(grepl("N", muts$mutation)), ]
muts$group = sapply(muts$mutation, FUN = function(x) switch(x, "AT" = 1, "TA" = 1, "AG" = 2, "TC" = 2, "AC" = 3, "TG" = 3, "GC" = 4, "CG" = 4, "GT" = 5, "CA" = 5, "GA" = 6, "CT" = 6))

for (group in unique(muts$group)) {
    current_base_df = data.frame(position = unique(total$position))
    print(group)
    for (base in unique(muts[muts$group == group, "ancestral"])) {
        c_base = paste("comptage_", base, sep = "")
        current_base_df[, c_base] = 0
        current_base_df[match(bases[bases$base == base, "position"], current_base_df$position), c_base] = bases[bases$base == base, "comptage"]
        for (bdest in unique(muts[muts$ancestral == base & muts$group == group, "reference"])) {
            mut = paste(base, bdest, sep = "")
            current_base_df[, mut] = 0
            current_base_df[match(muts[muts$mutation == mut, "position"], current_base_df$position), mut] = muts[muts$mutation == mut, "comptage"]
        }
    }
    if (!is.null(ncol(current_base_df[, grep("comptage", colnames(current_base_df))])) && !is.null(ncol(current_base_df[, grep("^[ACGT]+?$", colnames(current_base_df))]))) {
        total_both = rowSums(current_base_df[, grep("comptage", colnames(current_base_df))])
        current_base_df$relative_both = rowSums(current_base_df[, grep("^[ACGT]+?$", colnames(current_base_df))]) / rowSums(current_base_df[, grep("comptage", colnames(current_base_df))])
    } else {
        total_both = current_base_df[, grep("comptage", colnames(current_base_df))]
        current_base_df$relative_both = current_base_df[, grep("^[ACGT]+?$", colnames(current_base_df))] / current_base_df[, grep("comptage", colnames(current_base_df))]
    }
    current_base_df$error_both = (2 * sqrt(current_base_df$relative_both * total_both * (1 - current_base_df$relative_both))) / total_both
    current_base_df$mean10 = mean10pb(current_base_df$relative_both)
    current_base_df$error10 = mean10pb(current_base_df$error_both)
    write_tsv(current_base_df[current_base_df$position %in% -220:5000, ], paste("group", group, "-", mut, "_and_reverse.tsv", sep = ""))
}

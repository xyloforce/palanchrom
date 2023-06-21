library(readr)
library(ggplot2)
library(stringr)
library(cowplot)
library(extrafont)

args = commandArgs(trailingOnly=TRUE)
#args = c("palanchrom/data/2022_03_24_results/2022_03_24_hg38_CG_formatted/", "palanchrom/data/2022_03_24_results/2022_03_24_panTro5_CG_formatted/", "human_snps/results_20220405/formatted_CG/", "chimp_snps/results/02_02_22_panTro5_CG_formatted/", "SNPxNIEBnormal.png")

source("~/setThemePoster.R")
file_contents = c("no TE", "FP Alu")
min_win = -25
max_win = 120

# normalise_muts = function(df) {return(df)}
normalise_muts = function(df) {
    df = as.data.frame(df)
    for(value in unique(df$source)) {
        df[df$source == value,"mean10"] = df[df$source == value,"mean10"] / mean(df[df$source == value,"mean10"], na.rm = TRUE)
    }
    return(df)
}
#### PLOT SOURCE X DEST X CONTEXT ################################################################################################

## load muts according to wildcard

print("loading data from first folder")

filelist = list.files(path = args[1], "*_and_reverse.tsv", full.names = TRUE)

df = read_tsv(filelist[1], show_col_types = FALSE)
filename = tail(strsplit(filelist[1], split = "/")[[1]], n = 1)
df$source = strsplit(filename, split = "_|-")[[1]][2]
df = df[,c("position", "source", "mean10")]

count = 1
for(filepath in filelist[2:length(filelist)]) {
	tmp = read_tsv(filepath, show_col_types = FALSE)
	filename = tail(strsplit(filepath, split = "/")[[1]], n = 1)
    tmp$source = strsplit(filename, split = "_|-")[[1]][2]
    tmp = tmp[,c("position", "source", "mean10")]
	df = rbind(df, tmp)
}

# print(strsplit(filelist[1], split = "/")[[1]])
df$type = file_contents[1]

# df = normalise_muts(df)

## now we loop

for(arg in args[2:(length(args)-1)]) {
    print("loading data")

    filelist = list.files(path = arg, "*_and_reverse.tsv", full.names = TRUE)

    df2 = read_tsv(filelist[1], show_col_types = FALSE)
    filename = tail(strsplit(filelist[1], split = "/")[[1]], n = 1)
    df2$source = strsplit(filename, split = "_|-")[[1]][2]
    df2 = df2[,c("position", "source", "mean10")]

    for(filepath in filelist[2:length(filelist)]) {
        tmp = read_tsv(filepath, show_col_types = FALSE)
#         print(tail(strsplit(filepath, split = "/")[[1]], n = 1))
        filename = tail(strsplit(filepath, split = "/")[[1]], n = 1)
        tmp$source = strsplit(filename, split = "_|-")[[1]][2]
        tmp = tmp[,c("position", "source", "mean10")]
        df2 = rbind(df2, tmp)
    }
    count = count + 1
    df2$type = file_contents[count]
#     head(df2)
#    df2 = normalise_muts(df2)
    print("merging")

    df = rbind(df, df2)
}

write_tsv(df, "savestate.tsv")

print("plotting")

correct_labels <- c("TG" = "A→C - T→G (transversion)",
                    "TC" = "A→G - T→C (transition)",
                    "TA" = "A→T - T→A (transversion)",
                    "GT" = "C→A - G→T (transversion)",
                    "GC" = "C→G - G→C (transversion)",
                    "GA" = "C→T - G→A (transition)"
                    )

plot = ggplot(data = df, aes(x = position, y = mean10, color = type)) + facet_wrap( ~ source, scales = "free", labeller = as_labeller(correct_labels)) + geom_line(size = 1.5) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
	scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
	geom_vline(xintercept = c(0, 133), color = "black", size = 2) +
	xlab("Distance from the NIEB 3' border") +
	ylab("Mutation rate") +
	theme_poster +
	theme(strip.placement = "outside", legend.position = "bottom", legend.key.width = unit(1, "in")) +
	scale_color_manual(values = c("#009C7D", "#FF713D")) +
	coord_cartesian(xlim = c(min_win, max_win))

ggsave(args[length(args)], width = 20, height = 9)

library(ggplot2)
library(stringr)
library(cowplot)
library(extrafont)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")

# ARGS ARE :
# path1:path2(:path3...)
# type1:type2(:type3...)
# filename
# optionnal : min
# optionnal : max

file_paths = str_split(args[1], pattern = ":", simplify = TRUE)[1, ]
file_contents = str_split(args[2], pattern = ":", simplify = TRUE)[1, ]

min_win = -50
max_win = 300
if (length(args) > 3) {
    min_win = as.numeric(args[4])
    max_win = as.numeric(args[5])
}

normalise_muts = function(df) {
    for (value in unique(df$source)) {
        df[df$source == value, "mean10"] = df[df$source == value, "mean10"] /
            mean(df[df$source == value, "mean10"], na.rm = TRUE)
        df[df$source == value, "error10"] = df[df$source == value, "error10"] /
            mean(df[df$source == value, "error10"], na.rm = TRUE)
    }
    return(df)
}
####################### PLOT SOURCE X DEST X CONTEXT #######################

####################### load muts according to wildcard

print("loading data from first folder")

filelist = list.files(path = file_paths[1], "*_and_reverse.tsv",
                      full.names = TRUE)
filelist = c(paste(file_paths[1], "/total.tsv", sep = ""), filelist)
df = read.delim(filelist[1])
df$source = "Total"
df = df[, c("position", "source", "mean10", "error10")]

count = 1
for (filepath in filelist[2:length(filelist)]) {
    tmp = read.delim(filepath)
    filename = tail(strsplit(filepath, split = "/")[[1]], n = 1)
    tmp$source = strsplit(filename, split = "_|-")[[1]][2]
    tmp = tmp[, c("position", "source", "mean10", "error10")]
    df = rbind(df, tmp)
}

df$type = file_contents[1]

# df = normalise_muts(df)

## now we loop

for (arg in file_paths[2:(length(file_paths))]) {
    print("loading data")

    filelist = list.files(path = arg, "*_and_reverse.tsv", full.names = TRUE)
    filelist = c(paste(arg, "/total.tsv", sep = ""), filelist)
    df2 = read.delim(filelist[1])
    df2$source = "Total"
    df2 = df2[, c("position", "source", "mean10", "error10")]

    for (filepath in filelist[2:length(filelist)]) {
        tmp = read.delim(filepath)
        filename = tail(strsplit(filepath, split = "/")[[1]], n = 1)
        tmp$source = strsplit(filename, split = "_|-")[[1]][2]
        tmp = tmp[, c("position", "source", "mean10", "error10")]
        df2 = rbind(df2, tmp)
    }
    count = count + 1
    df2$type = file_contents[count]
    # df2 = normalise_muts(df2)
    print("merging")

    df = rbind(df, df2)
}

write.table(df, "savestate.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

print("plotting")

correct_labels <- c("TG" = "A→C - T→G (transversion)",
                    "TC" = "A→G - T→C (transition)",
                    "TA" = "A→T - T→A (transversion)",
                    "GT" = "C→A - G→T (transversion)",
                    "GC" = "C→G - G→C (transversion)",
                    "GA" = "C→T - G→A (transition)",
                    "Total" = "Global mutation rate"
                    )

vline_coords = c(0, 133)
vline_coords = vline_coords[vline_coords > min_win & vline_coords < max_win]

df[df$position %% 10 != 0, "error10"] = 0
df = df[df$position %in% min_win:max_win, ]
plot = ggplot(data = df, aes(x = position, y = mean10, color = type)) +
    facet_wrap(~ source, scales = "free",
               labeller = as_labeller(correct_labels)) +
    geom_line(linewidth = 1.5) +
    geom_vline(xintercept = vline_coords, color = "black", linewidth = 1) +
    geom_errorbar(aes(ymin = mean10 - error10, ymax = mean10 + error10),
                      color = "black") +
    xlab("position") +
    ylab("mutation rate") +
    theme_poster +
    scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
    scale_color_manual(values = c("#f53da4", "#f5a80c", "#0c65f5", "#22f518"))

ggsave(args[3], width = 20, height = 9)

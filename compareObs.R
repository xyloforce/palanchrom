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
normalize = FALSE
if (length(args) > 5) {
    normalize = TRUE
}

# delete_additionnal = FALSE
keep_groups = c()
if (length(args) > 6) {
    keep_groups = args[7]
    keep_groups = str_split(keep_groups, pattern = ":", simplify = TRUE)[1, ]
    keep_groups = c(keep_groups, "Total")
}

normalise_muts = function(df) {
    for (value in unique(df$source)) {
        mean_mut_rate = sum(df[df$source == value, "mutations"], na.rm = TRUE) /
                        sum(df[df$source == value, "comptage"], na.rm = TRUE)
        df[df$source == value, "mean10"] = df[df$source == value, "mean10"] /
            mean_mut_rate
        df[df$source == value, "error10"] = df[df$source == value, "error10"] /
            mean_mut_rate
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
df = df[, c("position", "comptage", "mutations",
            "source", "mean10", "error10")]

count = 1
for (filepath in filelist[2:length(filelist)]) {
    tmp = read.delim(filepath)
    filename = tail(strsplit(filepath, split = "/|\\\\")[[1]], n = 1)
    tmp$source = strsplit(filename, split = "-")[[1]][1]
    tmp$comptage = sum(tmp[, grep("comptage", colnames(tmp))])
    tmp$mutations = sum(tmp[, grep("[ACGT]{2}", colnames(tmp))])
    tmp = tmp[, c("position", "comptage", "mutations",
                  "source", "mean10", "error10")]
    df = rbind(df, tmp)
}

df$type = file_contents[1]

if (normalize) {
    df = normalise_muts(df)
}

## now we loop

for (arg in file_paths[2:(length(file_paths))]) {
    print("loading data")
    filelist = list.files(path = arg, "*_and_reverse.tsv", full.names = TRUE)
    filelist = c(paste(arg, "total.tsv", sep = ""), filelist)
    df2 = read.delim(filelist[1])
    df2$source = "Total"
    df2 = df2[, c("position", "comptage", "mutations",
                  "source", "mean10", "error10")]
    for (filepath in filelist[2:length(filelist)]) {
        tmp = read.delim(filepath)
        filename = tail(strsplit(filepath, split = "/|\\\\")[[1]], n = 1)
        tmp$source = strsplit(filename, split = "-")[[1]][1]
        tmp$comptage = sum(tmp[, grep("comptage", colnames(tmp))])
        tmp$mutations = sum(tmp[, grep("[ACGT]{2}", colnames(tmp))])
        tmp = tmp[, c("position", "comptage", "mutations",
                      "source", "mean10", "error10")]
        df2 = rbind(df2, tmp)
    }
    count = count + 1
    df2$type = file_contents[count]
    if (normalize) {
        df2 = normalise_muts(df2)
    }
    print("merging")
    df = rbind(df, df2)
}

df = df[df$position > min_win & df$position < max_win, ]
write.csv2(df, "savestate.csv", row.names = FALSE)

print("plotting")

correct_labels <- c("group1" = "A→C - T→G (transversion)",
                    "group4" = "A→G - T→C (transition)",
                    "group2" = "C→A - G→T (transversion)",
                    "group6" = "C→G - G→C (transversion)",
                    "group5" = "C→T - G→A (transition)",
                    "group3" = "A→T - T→A (transversion)",
                    "Total" = "Global mutation rate"
                    )

vline_coords = c(0, 133)
vline_coords = vline_coords[vline_coords > min_win & vline_coords < max_win]

spacing_bars = 30
for (x in seq_along(unique(df$type))) {
    type = unique(df$type)[x]
    df[(df$position +
        (spacing_bars / length(unique(df$type)) * x)) %%
        spacing_bars != 0 &
        df$type == type, "error10"] = NA
        # ensure that bars don't fall at the same place for different types
}

if (normalize) {
    label_axis = "relative mutation rate"
} else {
    label_axis = "mutation rate"
}

if(length(keep_groups) > 0) {
    head(df)
    df = df[df$source %in% keep_groups,]
    head(df)
}

# if (delete_additionnal) {
#     types = unique(df$type)
#     elements = unique(df[df$type == types[1], "source"])
#     print(elements)
#     for (type in types[2:length(types)]) {
#         elements = intersect(elements, unique(df[df$type == type, "source"]))
#     }
#     print(elements)
#     df = df[df$source %in% elements, ]
#     print(head(df))
# }

df[df$position < min_win + 10, "error10"] = NA
df[df$position > max_win - 10, "error10"] = NA
df = df[df$position %in% min_win:max_win, ]
plot1 = ggplot(data = df[df$source != "Total", ],
               aes(x = position, y = mean10, color = type)) +
    facet_wrap(~ source, scales = "free",
               labeller = as_labeller(correct_labels)) +
    geom_line(linewidth = 1.5) +
    geom_vline(xintercept = vline_coords, linewidth = 1,
               color = "black") +
    geom_errorbar(aes(ymin = mean10 - error10, ymax = mean10 + error10),
                  width = .2) +
    xlab("position") +
    ylab(label_axis) +
    theme_poster +
    scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
    # scale_color_manual(values = c("#f53da4", "#f5a80c", "#0c65f5", "#22f518")) +
    theme(strip.text = element_text(size = 19))
plot2 = ggplot(data = df[df$source == "Total", ],
               aes(x = position, y = mean10, color = type)) +
    facet_wrap(~ source, scales = "free",
               labeller = as_labeller(correct_labels)) +
    geom_line(linewidth = 1.5) +
    geom_vline(xintercept = vline_coords, color = "black", linewidth = 1) +
    geom_errorbar(aes(ymin = mean10 - error10, ymax = mean10 + error10),
                  width = .2) +
    xlab("position") +
    ylab(label_axis) +
    theme_poster +
    scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) 
    # scale_color_manual(values = c("#f53da4", "#f5a80c", "#0c65f5", "#22f518"))

figure = plot_grid(plot1, plot2, labels = c("A", "B"), rel_widths = c(0.6, 0.3))

ggsave(args[3], width = 24, height = 9)

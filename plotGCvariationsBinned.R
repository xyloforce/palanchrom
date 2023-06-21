library(zoo)
library(ggh4x)
library(stringr)
library(ggplot2)
library(extrafont)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")
# args = c("formattedNCG", "gc_mutations.png", 20)

bin_size = 20
reverse = FALSE
xlimM = -50
xlimP = 270
ylimP = 0.008
ylimM = 0

comptages = FALSE
vertical_bars = FALSE

print("loading df muts")
fileList = list.files(path = args[1], "*_mutations.tsv", full.names = TRUE)

df = data.frame()
for (filename in fileList) {
    print(filename)
    df2 = read.delim(filename)
    if (grepl("A_mutations|T_mutations", filename)) { # want mutations going in the GC direction
        df = rbind(df, cbind(df2$position, "W → S", str_match(filename, "([AT])_mutations.tsv")[, 2], df2$comptage, rowSums(df2[, grepl("^[AT][GC]$", colnames(df2))], na.rm = TRUE)))
        df = rbind(df, cbind(df2$position, "W → W | S → S", str_match(filename, "([AT])_mutations.tsv")[, 2], df2$comptage, df2[, grepl("^[AT][AT]$", colnames(df2))]))
    } else { # want mutations going in the AT direction
        df = rbind(df, cbind(df2$position, "S → W", str_match(filename, "([GC])_mutations.tsv")[, 2], df2$comptage, rowSums(df2[, grepl("^[GC][AT]$", colnames(df2))], na.rm = TRUE)))
        df = rbind(df, cbind(df2$position, "W → W | S → S", str_match(filename, "([GC])_mutations.tsv")[, 2], df2$comptage, df2[, grepl("^[GC][GC]$", colnames(df2))]))
    }
}
colnames(df) = c("position", "type", "base", "total", "count")

df = aggregate(df[, c("total", "count")], by = list(position = df$position, type = df$type), FUN = function(x) {
    return(sum(as.numeric(x), na.rm = TRUE))
})
df$position = as.numeric(df$position)
df = df[df$position < xlimP + bin_size & df$position > xlimM - bin_size, ]
df = df[order(df$position), ]

print("Binning")
df[df$type == "W → S", "total_bases"] = rollapply(df[df$type == "W → S", "total"], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")
df[df$type == "S → W", "total_bases"] = rollapply(df[df$type == "S → W", "total"], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")
df[df$type == "W → W | S → S", "total_bases"] = rollapply(df[df$type == "W → W | S → S", "total"], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")
df[df$type == "W → S", "slidedCount"] = rollapply(df[df$type == "W → S", "count"], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")
df[df$type == "S → W", "slidedCount"] = rollapply(df[df$type == "S → W", "count"], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")
df[df$type == "W → W | S → S", "slidedCount"] = rollapply(df[df$type == "W → W | S → S", "count"], width = bin_size, FUN = sum, na.rm = TRUE, fill = NA, align = "center")

print("Calculating ratios")
# totalMut = aggregate(df[,"count"], by = list(bin = df$bin), FUN = function(x) {return(sum(as.numeric(x), na.rm = TRUE))})
totalMut = aggregate(df[, "slidedCount"], by = list(position = df$position), FUN = function(x) {
    return(sum(as.numeric(x), na.rm = TRUE))
})

# df$ratioMut = df$count / totalMut[match(df$bin, totalMut$bin), "x"]
# df$ratioMut = df$slidedCount / totalMut[match(df$position, totalMut$position), "x"]
# df$ratio = df$count / df$total
df$ratio = df$slidedCount / df$total_bases

df$error = (sqrt((df$ratio * (1 - df$ratio)) / df$total))
df[is.na(df$error), "error"] = 1
df[df$pos %% 10 != 0, "error"] = NA
df$ymax = df$ratio + df$error
df$ymin = df$ratio - df$error

save(df, file = "df.Rdata")

if (substr(args[1], nchar(args[1]), nchar(args[1])) == "/") {
    args[1] = substr(args[1], 1, nchar(args[1]) - 1)
}

final_df = reshape(data = df, idvar = "position", varying = c("ratio", "total_bases"), v.name = "value", times = c("ratio", "total_bases"), new.row.names = 1:1000000, direction = "long")
head(final_df)

save(final_df, file = "final_df.Rdata")

color = c("W → S" = "#e66e00", "S → W" = "#17e629", "W → W | S → S" = "#4e17e6")
# fill = c("W → S" = "#ea9546", "S → W" = "#5eea69", "W → W | S → S" = "#835eea")

if(comptages) {
    plot = ggplot(data = final_df, aes(x = position, y = value, color = type)) + geom_line(linewidth = 1) + geom_errorbar(aes(ymin = ymin, ymax = ymax, color = type)) + theme_poster + theme(axis.text.x = element_text(angle = 90)) + scale_color_manual(values = color) + facet_grid(rows = vars(time), scales = "free") + xlab("position relative to the NIEB") + ylab(paste("rate of mutation smoothed over ", bin_size, " bp", sep = ""))
} else {
    final_df = final_df[final_df$time == "ratio",]
    plot = ggplot(data = final_df, aes(x = position, y = value, color = type)) + geom_line(linewidth = 1) + geom_errorbar(aes(ymin = ymin, ymax = ymax, color = type)) + theme_poster + theme(axis.text.x = element_text(angle = 90)) + scale_color_manual(values = color) + xlab("Distance from the barrier") + ylab(paste("Mutation rate over ", bin_size, " bp", sep = ""))
}

if(vertical_bars) {
	plot = plot + geom_vline(xintercept = c(0, 133, 266), color = "black") + geom_vline(xintercept = c(140), linetype = 2, color = "red")
}

if (reverse) {
    plot = plot + scale_x_reverse()
}

plot = plot + facetted_pos_scales(y = list(scale_y_continuous(limits = c(ylimM, ylimP)), NULL))

ggsave(args[2], width = 10, height = 9)

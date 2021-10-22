library(readr)
library(ggplot2)


countCPG = read_tsv("Scripts/data/long_countBasesCPG.tsv", col_names = FALSE)
countCPG$type = "CPG"
names(countCPG) = c("pos", "base", "count", "type")
subset = countCPG[countCPG$pos > -100 & countCPG$pos < 500,]
rm(countCPG)

countNCPG = read_tsv("Scripts/data/long_countBasesNCPG.tsv", col_names = FALSE)
countNCPG$type = "NCPG"
names(countNCPG) = c("pos", "base", "count", "type")
subset = rbind(subset, countNCPG[countNCPG$pos > -100 & countNCPG$pos < 500,])
rm(countNCPG)

ancestral = sapply(subset$base, function(x) {
  if(x %in% c("A", "B", "D", "E", "F")) {
    return("A")
  } else if (x %in% c("H", "C", "I", "J", "K")) {
    return("C")
  } else if(x %in% c("L", "M", "G", "O", "P")) {
    return("G")
  } else if(x %in% c("Q", "R", "S", "T", "V")) {
    return("T")
  }
})

subset = cbind(subset, ancestral)
rm(ancestral)

# check if there is a trend in mutation around barriers
sumByPos = aggregate(subset$count, by=list(subset$pos), FUN=sum)
percentByPos = apply(subset, MARGIN = 1, function(x) {
  return (as.numeric(x[3]) / sumByPos[sumByPos$`Group.1` == as.numeric(x[1]),"x"] * 100)
})
subset = cbind(subset, percentByPos)
isMut = sapply(subset$base, function(x) {
  if(x %in% c("A", "C", "G", "T", "K")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
subset = cbind(subset, isMut)

muted = aggregate(subset$count, by=list(subset$pos, subset$isMut), FUN = sum)
per = apply(muted, MARGIN = 1, function(x) {
  return (x[3]/sum(muted[muted$Group.1 == x[1],"x"])*100)
})
muted = cbind(muted, per)
names(muted) = c("pos", "isMut", "count", "percent")
ggplot(muted, aes(x=pos, y=percent, fill=isMut)) + geom_bar(stat="identity")

# get mutations rates by cells
ggplot(subset, aes(x=pos, y=count, color = type)) + geom_point() + facet_wrap(~base) + scale_y_log10()

# check if we can get back the trend in GC around barriers
subsetCPG = subset[subset$type=="CPG",]
values = aggregate(subsetCPG$count, by=list(subsetCPG$pos), FUN = sum)
ggplot(values, aes(x=Group.1, y = x)) + geom_point()

## and in percent ?
subset = data[data$pos < 500 & data$pos > -100,]
count = aggregate(subset$count, by=list(subset$pos, subset$type), FUN = sum)
percentCPG = apply(count, MARGIN = 1, function(x) {
  return(as.numeric(x[3])/sum(count[as.numeric(x[1]) == count$Group.1,"x"])*100)
})
count$percentCPG = percentCPG
names(count) = c("pos", "type", "count", "percent")
ggplot(count[count$type == "CPG",], aes(x=pos, y = percent)) + geom_smooth()
ggplot(count[count$type == "CPG",], aes(x=pos, y = percent)) + geom_point()

################################################
data = read_tsv("Scripts/data/long_countBasesNCPG.tsv", col_names = FALSE)
data$type = "NCPG"
names(data) = c("pos", "base", "count", "type")
data = rbind(data, countCPG)
rm(countCPG)

ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count, color = type)) + geom_point(alpha =0.5) + facet_wrap(~base)
ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count, color = base)) + geom_point() + facet_wrap(~type) + scale_y_log10()
ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count, shape = base)) + geom_point() + facet_wrap(~type) + scale_y_log10()
ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count)) + geom_point() + facet_wrap(base~type) + scale_y_log10()
ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count, shape = ancestral, color = base)) + geom_point() + facet_wrap(~type) + scale_y_log10()
ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count, color = base, shape = type)) + geom_point() + facet_wrap(~ancestral) + scale_y_log10()

ggplot(data[data$pos > -50 & data$pos < 50,], aes(x= pos, y=count, fill = type)) + geom_bar(stat = "identity") + scale_y_log10()

subset = data[data$pos > -50 & data$pos < 50,]
sumByPos = aggregate(subset$count, by=list(subset$pos), FUN=sum)
percentByPos = apply(subset, MARGIN = 1, function(x) {
  return (as.numeric(x[3]) / sumByPos[sumByPos$`Group.1` == as.numeric(x[1]),"x"] * 100)
})
subset = cbind(subset, percentByPos)
isMut = sapply(subset$base, function(x) {
  if(x %in% c("A", "C", "G", "T", "K")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
subset = cbind(subset, isMut)

muted = aggregate(subset$count, by=list(subset$pos, subset$isMut), FUN = sum)
per = apply(muted, MARGIN = 1, function(x) {
  return (x[3]/sum(muted[muted$Group.1 == x[1],"x"])*100)
})
muted = cbind(muted, per)
names(muted) = c("pos", "isMut", "count", "percent")
ggplot(muted, aes(x=pos, y=percent, color=isMut)) + geom_point()

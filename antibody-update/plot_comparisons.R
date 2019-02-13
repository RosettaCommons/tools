# plot antibody database comparisons
library(ggplot2)
library(reshape2)

# hard coded dir
setwd("~/Rosetta/tools/antibody-update/")

cdr.names <- c("h1", "h2", "h3", "l1", "l2", "l3")

# read sequence comparisons
temp.df <- data.frame()
loop.sequences <- data.frame()
for (name in cdr.names) {
  temp.df <- read.csv(paste0("comparison/", name, "_sequences.csv"), header=T, stringsAsFactors = F)
  loop.sequences <- rbind(loop.sequences, temp.df)
}

# read length comparisons
temp.df <- data.frame()
loop.lengths <- data.frame()
for (name in cdr.names) {
  temp.df <- read.csv(paste0("comparison/", name, "_lengths.csv"), header=T, stringsAsFactors = F)
  loop.lengths <- rbind(loop.lengths, temp.df)
}
# some dumb data manipulation because I am not smart
loop.lengths <- melt(loop.lengths, id.vars = c("cdr", "length"))

# read B factors
temp.df <- data.frame()
bfactors <- data.frame()
for (name in cdr.names) {
  temp.df <- read.csv(paste0("comparison/", name, "_bfactors.csv"), header=T, stringsAsFactors = F)
  bfactors <- rbind(bfactors, temp.df)
}
# some dumb data manipulation because I am not smart
bfactors <- melt(bfactors, id.vars = c("cdr", "pdb"))

# read OCD comparisons -- 60 MB file, slow read
ocds <- read.csv("comparison/ocds.csv", header=T)
# and if it wasn't big enough, lets make it bigger!
ocds$diff <- ocds$old_ocd - ocds$new_ocd

# plot lengths -- exclude super long H3s (save as 8x6")
lengths.plot <- ggplot(subset(loop.lengths,length<25), aes(x=length,y=value,fill=variable)) + geom_bar(stat="identity", position="dodge") + facet_wrap(~cdr, scales="free") + guides(fill=F) + theme_bw(base_size = 12) + xlab("Length") + ylab("Count") + ggtitle("CDR Length Distributions New (Blue) vs. Old (Red)")

# plot count of sequence mismatches 
sequences.plot <- ggplot(loop.sequences,aes(x=cdr)) + geom_bar(fill="gray50",color="black") + theme_bw(base_size = 12) + xlab("CDR") + ylab("Count") + ggtitle("CDR Sequence Mismatches")

# plot Bfactors
bfactors.plot <- ggplot(bfactors, aes(x=value, fill=variable)) + geom_bar(position="dodge", color="black") + facet_wrap(~cdr, scales="free") + guides(fill=F) + theme_bw(base_size = 12) + xlab("B-Factor Value") + ylab("Count") + ggtitle("CDR B-Factor Values New (Blue) vs. Old (Red)")

# plot OCDs direct comparison
ocd.line.plot<-ggplot(subset(ocds, abs(diff)>4),aes(x=old_ocd,y=new_ocd)) + geom_point() + geom_abline(slope=1, intercept=0) + theme_bw(base_size = 12) + xlab("Old OCD") + ylab("New OCD") + ggtitle("OCD Comparison of Antibody Pairs")
# save due to size
#ggsave(filename = "~/Desktop/OCD_line.pdf", plot = ocd.line.plot, width = 6, height = 4, units = "in")

# plot density of OCD deltas, exclude around 0 because there are lots
ocd.delta.plot <- ggplot(subset(ocds, abs(diff)>4),aes(x=diff)) + geom_histogram(bins=100, fill="gray50", color="black") + theme_bw(base_size = 12) + xlab("Old - New OCD") + ylab("Density") + ggtitle("OCD Comparison of Antibody Pairs") + scale_y_log10(breaks=c(1,10,100,1000,10000))

# bonus plot OCD distributions
# read new angle metrics
angles.new <- read.csv("~/Rosetta/tools/antibody-update/info/angles.info", header=T)
angles.new$label <- "new"

# read old metrics
angles.old <- read.csv("~/Rosetta/tools/antibody/angles.sc", header=T)
angles.old$label <- "old"
names(angles.old) <- c("pdb", "distance", "heavy_opening_angle", "light_opening_angle", "packing_angle", "label")

# plot each angle metric distribution
# these are intended to be saved as 6 by 4 inch plots!
angles <- rbind(angles.new, angles.old)

# distance
dist.plot <- ggplot(angles, aes(x=heavy_opening_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + xlim(13,17)  + guides(color=F) + ggtitle("Heavy Opening Angle in Old (Blue) vs. New (Red) AbDb") + xlab("H-Opening Angle (Deg)") + ylab("Density")

# heavy angle
hopen.plot <- ggplot(angles, aes(x=heavy_opening_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + guides(color=F) + ggtitle("Heavy Opening Angle in Old (Blue) vs. New (Red) AbDb") + xlab("H-Opening Angle (Deg)") + ylab("Density")

# light angle
lopen.plot <- ggplot(angles, aes(x=light_opening_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + guides(color=F) + ggtitle("Light Opening Angle in Old (Blue) vs. New (Red) AbDb") + xlab("L-Opening Angle (Deg)") + ylab("Density")

# packing angle
pa.plot <- ggplot(angles, aes(x=packing_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + guides(color=F) + ggtitle("Packing Angle in Old (Blue) vs. New (Red) AbDb") + xlab("Packing Angle (Deg)") + ylab("Density")

# print averages/sds
#print(sapply(angles.new, mean))
#print(sapply(angles.new, sd))

#print(sapply(angles.old, mean))
#print(sapply(angles.old, sd))

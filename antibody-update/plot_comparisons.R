# plot antibody database comparisons
library(ggplot2)
library(reshape2)

# hard coded dir -- set to your own!
setwd("~/Rosetta/tools/antibody-update/")

# if you want graft analysis, set to true
graft.analysis <- T
# if you want info file comparison, set to true
info.analysis <- T

cdr.names <- c("h1", "h2", "h3", "l1", "l2", "l3")

if (graft.analysis) {
  # read graft RMSDs
  new <- read.csv("new_db_grafts.csv", header=T)
  new$label <- "new"
  
  old <- read.csv("old_db_grafts.csv", header=T)
  old$label <- "old"
  
  # manipulate
  both <- rbind(old,new)
  melt.both <-melt(both, id.vars=c("pdb","model","label"))
  merge.both <- merge(subset(melt.both, label=="old"), subset(melt.both, label=="new"), by=c("pdb","model","variable"), all = T)
  
  # plot
  violin.plot<-ggplot(subset(mboth, model==0 & variable!="ocd"), aes(x=variable, y=value, fill=label)) + geom_violin() + theme_bw(base_size=12) + guides(fill=F) + ylab("RMSD (A)") + xlab("Region")
  ggsave(filename = "region_vs_rmsd_violins.pdf", plot = violin.plot,  width = 6, height = 4, units = "in")
  
  density.plot<-ggplot(subset(merge.both, model==0 ), aes(x=value.x - value.y)) + geom_density() + facet_wrap(~variable, scales="free") + theme_bw(base_size=12) + ylab("Density") + xlab("Old-New")
  ggsave(filename = "region_vs_rmsd_delta_density.pdf", plot = density.plot,  width = 6, height = 6, units = "in")
  
  paired.plot<-ggplot(subset(merge.both, model==0 ), aes(x=value.x, y=value.y)) + geom_point() + geom_abline() + geom_hline(yintercept = 1) + geom_vline(xintercept = 1) + facet_wrap(~variable, scales="free") + theme_bw(base_size=12) + ylab("New Value") + xlab("Old Value")
  ggsave(filename = "paired_comparison_model_0.pdf", plot = paired.plot,  width = 6, height = 6, units = "in")
}

# lots of hard code in here
# TODO: fix
if (info.analysis) {
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
  ggsave(filename = "database_lengths.pdf", plot = lengths.plot,  width = 8, height = 6, units = "in")
  
  # plot count of sequence mismatches 
  sequences.plot <- ggplot(loop.sequences,aes(x=cdr)) + geom_bar(fill="gray50",color="black") + theme_bw(base_size = 12) + xlab("CDR") + ylab("Count") + ggtitle("CDR Sequence Mismatches")
  ggsave(filename = "sequence_mismatch_counts.pdf", plot = sequences.plot,  width = 6, height = 4, units = "in")
  
  # plot Bfactors
  # confusing plot -- shows only the bfactor changes
  bfactors.plot <- ggplot(bfactors, aes(x=value, fill=variable)) + geom_bar(position="dodge", color="black") + facet_wrap(~cdr, scales="free") + guides(fill=F) + theme_bw(base_size = 12) + xlab("B-Factor Value") + ylab("Count") + ggtitle("CDR B-Factor Values New (Blue) vs. Old (Red)")
  ggsave(filename = "bfactor_change_counts.pdf", plot = bfactors.plot,  width = 6, height = 4, units = "in")
  
  # plot OCDs direct comparison
  ocd.line.plot<-ggplot(subset(ocds, abs(diff)>4),aes(x=old_ocd,y=new_ocd)) + geom_point() + geom_abline(slope=1, intercept=0) + theme_bw(base_size = 12) + xlab("Old OCD") + ylab("New OCD") + ggtitle("OCD Comparison of Antibody Pairs")
  # save due to size
  ggsave(filename = "OCD_line.pdf", plot = ocd.line.plot, width = 6, height = 4, units = "in")
  
  # plot density of OCD deltas, exclude around 0 because there are lots
  #ocd.delta.plot <- ggplot(subset(ocds, abs(diff)>4),aes(x=diff)) + geom_histogram(bins=100, fill="gray50", color="black") + theme_bw(base_size = 12) + xlab("Old - New OCD") + ylab("Density") + ggtitle("OCD Comparison of Antibody Pairs") + scale_y_log10(breaks=c(1,10,100,1000,10000))
  
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
  dist.plot <- ggplot(angles, aes(x=distance, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + xlim(13,17)  + guides(color=F) + ggtitle("VH-VL Distances in Old (Blue) vs. New (Red) AbDb") + xlab("Distance (A)") + ylab("Density")
  ggsave(filename = "OCD_distance_dist.pdf", plot = dist.plot, width = 6, height = 4, units = "in")
  
  # heavy angle
  hopen.plot <- ggplot(angles, aes(x=heavy_opening_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + guides(color=F) + ggtitle("Heavy Opening Angle in Old (Blue) vs. New (Red) AbDb") + xlab("H-Opening Angle (Deg)") + ylab("Density")
  ggsave(filename = "OCD_HOA_dist.pdf", plot = hopen.plot, width = 6, height = 4, units = "in")
  
  # light angle
  lopen.plot <- ggplot(angles, aes(x=light_opening_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + guides(color=F) + ggtitle("Light Opening Angle in Old (Blue) vs. New (Red) AbDb") + xlab("L-Opening Angle (Deg)") + ylab("Density")
  ggsave(filename = "OCD_LOA_dist.pdf", plot = lopen.plot, width = 6, height = 4, units = "in")
  
  # packing angle
  pa.plot <- ggplot(angles, aes(x=packing_angle, color=label)) + geom_density(size=1) + theme_bw(base_size=12) + guides(color=F) + ggtitle("Packing Angle in Old (Blue) vs. New (Red) AbDb") + xlab("Packing Angle (Deg)") + ylab("Density")
  ggsave(filename = "OCD_PA_dist.pdf", plot = pa.plot, width = 6, height = 4, units = "in")
  
  # print averages/sds
  print("New")
  print(sapply(angles.new[,-c(1,6)], mean))
  print(sapply(angles.new[,-c(1,6)], sd))
  
  print("Old")
  print(sapply(angles.old[,-c(1,6)], mean))
  print(sapply(angles.old[,-c(1,6)], sd))  
}


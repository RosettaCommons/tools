#!/usr/bin/env Rscript
# kink_plots.R
# plot geometry of the H3 kink for a set of modeled antibody structures
# J. Gray June 2013


nd12=read.table("ND2-12.H3flex",header=T)
nd9=read.table("ND2-9.H3flex",header=T)
nd15=read.table("ND2-15.H3flex",header=T)
jqx=read.table("91215.hom.bb.H3flex",header=T)

library(ggplot2)
p = qplot(jqx$H3len,jqx$H3rmsd,data=jqx)
p + geom_abline()

nd2 = rbind(nd9,nd12,nd15)
names(nd2)

nd2$set = factor("ND2")
jqx$set = factor("PDB")
a=as.character(nd2$Ab)
nd2$Ab = a
dat=rbind(jqx,nd2)

p = ggplot(dat,aes(factor(set,H3len),rmsd))
names(dat)[names(dat)=="H3rmsd"] <- "H3flex"
p + geom_boxplot(aes(fill=factor(set))) + theme_grey(base_size=20) + xlab("CDR H3 loop length") + ylab("Mean H3 RMSD in Top10 Models (Angstrom)")


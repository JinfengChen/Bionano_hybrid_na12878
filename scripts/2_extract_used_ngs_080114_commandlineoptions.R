#!/usr/bin/env Rscript

source("/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/scripts/fileIO.R")
args <- commandArgs(TRUE)

     #R --no-save -args x y z < script.R


all.cmap.m <- read.cmap(args[1])

ids <- read.table(args[2], header = T, sep="\t")

out.file.name <- args[3]

selected.ids <- ids[ids$newCmapId != -1,]
selected.cmap.m = all.cmap.m[all.cmap.m[, "CMapId"] %in% selected.ids$newCmapId,]

write.cmap(cmap.matrix = selected.cmap.m, file = out.file.name)


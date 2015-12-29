#!/usr/bin/env Rscript

source("fileIO.R")
args <- commandArgs(TRUE)

#from oscar's R scripts, useful execution conventions
#R --no-save -args x y z < script.R

R --no-save --args /sc/orga/scratch/pendlm02/testBioNano/end2end/end_to_end/new_quivered_assembly/consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels.cmap \
/sc/orga/scratch/pendlm02/testBioNano/end2end/end_to_end/align_usedQuiveredCelera_v2Hybrid/used.translation.fasta.cmap.txt \
/sc/orga/scratch/pendlm02/testBioNano/end2end/end_to_end/align_usedQuiveredCelera_v2Hybrid/used.newSeq.cmap <  2_extract_used_ngs_080114_commandlineoptions.R

#all.cmap.m <- read.cmap("/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/new_pac_bio_seq/consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels.cmap") 

all.cmap.m <- read.cmap(args[1])

#ids <- read.table("/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/used.translation.fasta.cmap.txt", header = T)
ids <- read.table(args[2], header = T, sep="\t")

#out.file.name <- "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/used.newSeq.cmap";
out.file.name <- args[3]

selected.ids <- ids[ids$newCmapId != -1,]
selected.cmap.m = all.cmap.m[all.cmap.m[, "CMapId"] %in% selected.ids$newCmapId,]
#selected.cmap.m = all.cmap.m[colnames(all.cmap.m)$CMapId  %in% selected.ids$newCmapId,]
write.cmap(cmap.matrix = selected.cmap.m, file = out.file.name)


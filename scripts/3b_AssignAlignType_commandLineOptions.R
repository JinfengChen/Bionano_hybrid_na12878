source("fileIO.R")

args <- commandArgs(TRUE)

# input, should be the alignment and the two cmap files used in the alignment
#xmap.m = read.xmap("/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid.xmap")
xmap.m = read.xmap(args[1])

#ngs.cmap.m = read.cmap("/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_r.cmap")
ngs.cmap.m = read.cmap(args[2])

#bn.cmap.m = read.cmap("/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_q.cmap")
bn.cmap.m = read.cmap(args[3])
# output
#sticky.xmap.name = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/selectedNewV1Seq.vs.v2Hybrid.xmap.txt";
sticky.xmap.name = args[4]
#filtered.ngs.cmap.name = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/hybrid.filtered.cmap"
filtered.ngs.cmap.name = args[5]
#filtered.bionano.cmap.name = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/ngs.filtered.cmap"
filtered.bionano.cmap.name = args[6]
#overhang.count.name = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/sticky.xmap.overhang.count.txt"
overhang.count.name = args[7]

# Assume Ref is ngs
i = 1  # 164

max.overhang.lab = 5			# very stringent to deal with all contig, previously tried 3

lab.cnt.m = c()
for (i in 1:NROW(xmap.m)) {
    # we first find out the match
    this.xmap.row = xmap.m[i,]

    this.ngs.cmap.rows = ngs.cmap.m[ngs.cmap.m[,"CMapId"] == as.character(this.xmap.row["RefContigID"]),]
    this.bn.cmap.rows = bn.cmap.m[bn.cmap.m[,"CMapId"] == as.character(this.xmap.row["QryContigID"]),]
    

    # we shift 20bp to find out only labels extending outside alignment
    this.ngs.cmap.rows = this.ngs.cmap.rows[this.ngs.cmap.rows[, "LabelChannel"] == 1,]
    ref.left.lab.cnt = sum(as.numeric(this.ngs.cmap.rows[, "Position"]) < (as.numeric(this.xmap.row[, "RefStartPos"]) - 20) ) 
    ref.right.lab.cnt = sum(as.numeric(this.ngs.cmap.rows[, "Position"]) > (as.numeric(this.xmap.row[, "RefEndPos"]) + 20) ) 

    # Ref is always + while qry can be either +/-. In this case we use ngs.left.lab.cnt to annotate ref 
    # bn.left|right.lab.cnt to annotate number of labels left, right based on alignment orientation, aka reference orientation instead of query orientation. 
    this.bn.cmap.rows = this.bn.cmap.rows[this.bn.cmap.rows[, "LabelChannel"] == 1,]
    if (this.xmap.row["Orientation"] == "+") {
        qry.left.lab.cnt = sum(as.numeric(this.bn.cmap.rows[, "Position"]) < (as.numeric(this.xmap.row[, "QryStartPos"]) - 20) ) 
        qry.right.lab.cnt = sum(as.numeric(this.bn.cmap.rows[, "Position"]) > (as.numeric(this.xmap.row[, "QryEndPos"]) + 20) ) 
    } else {
        qry.right.lab.cnt = sum(as.numeric(this.bn.cmap.rows[, "Position"]) < (as.numeric(this.xmap.row[, "QryEndPos"]) - 20) ) 
        qry.left.lab.cnt = sum(as.numeric(this.bn.cmap.rows[, "Position"]) > (as.numeric(this.xmap.row[, "QryStartPos"]) + 20) )
    }
    
    # we use extend.max.lab to decide the max number of label can extend outside as outlier and still classify as a nested "subset" alignment
    # we use outlier.min.lab to decide if one side of alignment is enough to defined as an outlier.
    lab.cnt.m = rbind(lab.cnt.m, c("XmapEntryID" = as.character(this.xmap.row["XmapEntryID"]), 
        "ref.left.lab.cnt" = ref.left.lab.cnt, "ref.right.lab.cnt" = ref.right.lab.cnt, 
        "qry.left.lab.cnt" = qry.left.lab.cnt, "qry.right.lab.cnt" = qry.right.lab.cnt))
} 


T.cutoff = 15				# very stringent
lab.cnt.m.T =  lab.cnt.m[as.numeric(xmap.m[, "Confidence"]) >= T.cutoff,]

sticky.xmapIds = as.numeric(lab.cnt.m.T[((as.numeric(lab.cnt.m.T[, "ref.left.lab.cnt"]) > max.overhang.lab) & (as.numeric(lab.cnt.m.T[, "qry.left.lab.cnt"]) > max.overhang.lab)) | 
    ((as.numeric(lab.cnt.m.T[, "ref.right.lab.cnt"]) > max.overhang.lab) & (as.numeric(lab.cnt.m.T[, "qry.right.lab.cnt"]) > max.overhang.lab)), "XmapEntryID"])
sticky.lab.cnt.m.T = lab.cnt.m.T[lab.cnt.m.T[, "XmapEntryID"] %in% sticky.xmapIds,]
sticky.xmapIds
sticky.xmap = xmap.m[xmap.m[, "XmapEntryID"] %in% sticky.xmapIds,]
write.table(sticky.xmap, file=sticky.xmap.name, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(sticky.lab.cnt.m.T, file = overhang.count.name, sep="\t")

sticky.NGS.contigs = xmap.m[sticky.xmapIds, "RefContigID"]
sticky.NGS.contigs
sticky.BN.contigs = xmap.m[sticky.xmapIds, "QryContigID"]
sticky.BN.contigs

filtered.ngs.cmap = ngs.cmap.m[!ngs.cmap.m[, "CMapId"] %in% sticky.NGS.contigs,]
write.cmap(cmap.matrix = filtered.ngs.cmap, file = filtered.ngs.cmap.name)

filtered.bionano.cmap = bn.cmap.m[!bn.cmap.m[, "CMapId"] %in% sticky.BN.contigs,]
write.cmap(cmap.matrix = filtered.bionano.cmap, file = filtered.bionano.cmap.name)


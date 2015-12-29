# Nucmer is a dependency
SCRIPT_PATH="./"

cp $SCRIPT_PATH/scripts/fileIO.R .


#Step 0 (make the key files) 
#perl $SCRIPT_PATH/scripts//fa2cmap.pl -i input_data/consensus_diploid_for_andy_forncbi_filtered2b.fa  -n BspQI

#Step 1
perl $SCRIPT_PATH/scripts//1_find_usedCeleraId_in_newQuiveredCeleraId_031615_commandLineOptions.pl \
v1/output_results/mergeNGS_BN/ \
input_data/ \
mar3_NA12878.scf_BspQI_0Kb_0labels_key.txt \
./ \
input_data \
consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels_key.txt

#Step 2
R --no-save --args \
input_data/consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels.cmap \
align_usedQuiveredCelera_v2Hybrid/used.translation.fasta.cmap.txt \
align_usedQuiveredCelera_v2Hybrid/used.newSeq.cmap \
< $SCRIPT_PATH/scripts//2_extract_used_ngs_080114_commandlineoptions.R 

#Step 3 -woah! Be sure to edit the original xml files in /sc/orga/scratch/pendlm02/testBioNano/script_release/xml/ so that the
#              processor and memory count matches the available processors/memory
perl $SCRIPT_PATH/scripts//3_ngs_to_hybridscaffold_commandlineoptions.pl \
align_usedQuiveredCelera_v2Hybrid/ \
v2/output_results/mergeNGS_BN/ \
$SCRIPT_PATH/RefAligner2


mkdir align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/

#Step 3b
R --no-save --args \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid.xmap \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_r.cmap \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_q.cmap \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/selectedNewV1Seq.vs.v2Hybrid.xmap.txt \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/hybrid.filtered.cmap \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/ngs.filtered.cmap \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/sticky.xmap.overhang.count.txt \
< $SCRIPT_PATH/scripts/3b_AssignAlignType_commandLineOptions.R

#grep out bad xmap entries
cut -f 2-4 align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/selectedNewV1Seq.vs.v2Hybrid.xmap.txt > excluded_xmap_query_contigs.txt
grep -v -f excluded_xmap_query_contigs.txt \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid.xmap \
> reduced_selectedNewV1Seq.vs.v2Hybrid.xmap

# get overlapping contigs within a scaffold for final cleanup
perl scripts/identify_overlap_alignments_commandline.pl \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/ \
./ reduced_selectedNewV1Seq.vs.v2Hybrid.xmap 


# final step
# build fasta sequence
python $SCRIPT_PATH/scripts//findOverlapsBetweenContigs.py \
input_data/consensus_diploid_for_andy_forncbi_filtered2b.fa \
input_data/consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels_key.txt \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/fasta/ngs_overlap_status/adjacent_ngs_overlap_status.txt \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_r.cmap \
-a reduced_selectedNewV1Seq.vs.v2Hybrid.xmap  -p build_overlaps 2>&1 > /dev/null

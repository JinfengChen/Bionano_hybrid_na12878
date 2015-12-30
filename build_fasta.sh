#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=130gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`
#install lib using perl5.10.1
#../../perl5.10.1/perl-5.10.1/perl Build.PL --install_base /rhome/cjinfeng/BigData/software/Perl_lib
#./Build & ./Build install
#../../perl5.10.1/perl-5.10.1/perl Makefile.PL PREFIX=/rhome/cjinfeng/BigData/software/Perl_lib/
#make & make install
export PERL5LIB=/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/lib/:/rhome/cjinfeng/BigData/software/perl:/rhome/cjinfeng/BigData/software/Perl_lib/:/rhome/cjinfeng/BigData/software/Perl_lib/lib/site_perl/5.10.1/x86_64-linux/:/rhome/cjinfeng/BigData/software/Perl_lib/lib/site_perl/5.10.1:/rhome/cjinfeng/BigData/software/Perl_lib/lib/5.10.1:/rhome/cjinfeng/BigData/software/Perl_lib/lib/5.10.1/:/rhome/cjinfeng/BigData/software/Perl_lib/lib/perl5/x86_64-linux:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/scripts/perl5

#step1
#/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks scripts/1_find_usedCeleraId_in_newQuiveredCeleraId_031615_commandLineOptions.pl \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/v1/output_results/mergeNGS_BN \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/input_data \
#asm.scf_BspQI_0Kb_0labels_key.txt \
#./ \
#new_quivered_assembly \
#asm.scf_BspQI_0Kb_0labels_key.txt

#step2
#Rscript scripts/2_extract_used_ngs_080114_commandlineoptions.R --args \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/input_data/asm.scf_BspQI_0Kb_0labels.cmap \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/align_usedQuiveredCelera_v2Hybrid/used.translation.fasta.cmap.txt \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/align_usedQuiveredCelera_v2Hybrid/used.newSeq.cmap

#step3
#/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks scripts/3_ngs_to_hybridscaffold_commandlineoptions.pl \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/align_usedQuiveredCelera_v2Hybrid/ \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/v2/output_results/mergeNGS_BN/ \
#/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/RefAligner2
#mkdir /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/na12878_architecture/bionano/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/

#step3b
#R --no-save --args \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid.xmap \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_r.cmap \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_q.cmap \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/selectedNewV1Seq.vs.v2Hybrid.xmap.txt \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/hybrid.filtered.cmap \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/ngs.filtered.cmap \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/sticky.xmap.overhang.count.txt \
#< scripts/3b_AssignAlignType_commandLineOptions.R

#step4
#grep out bad xmap entries
#cut -f 2-4 align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/AssignAlignType/selectedNewV1Seq.vs.v2Hybrid.xmap.txt > excluded_xmap_query_contigs.txt
#grep -v -f excluded_xmap_query_contigs.txt \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid.xmap \
#> reduced_selectedNewV1Seq.vs.v2Hybrid.xmap

# get overlapping contigs within a scaffold for final cleanup
#/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks scripts/identify_overlap_alignments_commandline.pl \
#align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/ \
#./ reduced_selectedNewV1Seq.vs.v2Hybrid.xmap 

#step5
# build fasta sequence
module load perl
python scripts/findOverlapsBetweenContigs.py \
input_data/asm.scf.fasta \
input_data/asm.scf_BspQI_0Kb_0labels_key.txt \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/fasta/ngs_overlap_status/adjacent_ngs_overlap_status.txt \
align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/selectedNewV1Seq.vs.v2Hybrid_r.cmap \
/opt/linux/centos/7.x/x86_64/pkgs/mummer/3.23/bin \
-a reduced_selectedNewV1Seq.vs.v2Hybrid.xmap -p build_overlaps 2>&1 > /dev/null

#/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks rerun_coassembly_pipeline_citrus.pl
#/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks rerun_coassembly_pipeline_021615.pl
#~/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks  hybridScaffold_v1.pl -n input_data/Cclementina_v1.0_scaffolds.fa_contig.fasta -b input_data/Cclementina_v1.0_scaffolds_BspQI.cmap -c xml/hybridScaffold_config_v1.xml -o v1 -f
#~/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks  hybridScaffold_v1.pl -n input_data/NC_010473_mock_scaffolds.fna -b input_data/NC_010473_mock_scaffolds_BspQI.cmap -c xml/hybridScaffold_config_v1.xml -o v1 -f

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"


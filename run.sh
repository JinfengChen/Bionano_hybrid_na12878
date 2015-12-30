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
/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks rerun_coassembly_pipeline_citrus.pl
#/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks rerun_coassembly_pipeline_021615.pl
#~/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks  hybridScaffold_v1.pl -n input_data/Cclementina_v1.0_scaffolds.fa_contig.fasta -b input_data/Cclementina_v1.0_scaffolds_BspQI.cmap -c xml/hybridScaffold_config_v1.xml -o v1 -f
#~/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks  hybridScaffold_v1.pl -n input_data/NC_010473_mock_scaffolds.fna -b input_data/NC_010473_mock_scaffolds_BspQI.cmap -c xml/hybridScaffold_config_v1.xml -o v1 -f

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"


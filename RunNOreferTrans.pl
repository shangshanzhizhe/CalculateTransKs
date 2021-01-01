#!/usr/bin/env perl
#######Created by Duxin in LZU########
#######NO Refer Transcript 1.0########
use warnings;
use strict;

my ($trans,$anno) = @ARGV;;
die "Usage:\nperl $0 <Trininty assembly results> <annotation: yes or no(default)>\n" if(@ARGV<1);
$anno = "no" unless($anno);
$trans =~ s/^\S+(Trinity\.fasta)/$1/;
my @j = <./Trinity.fasta>;
`ln -s $trans` if(@j != 1);
my $uniprot_sprot_db="/data/00/user/user105/Trinotate_db/uniprot_sprot.pep";
my $pfamdb = "/data/00/user/user105/Trinotate_db";

print "###########01 delete splice and LongOrfs##########\n";
print "Deleting splice\n";
`perl ~/miniconda3/envs/Trinity/opt/trinity-2.1.1/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $trans > unigene.fasta 2</dev/null`;
`cd-hit-est -i unigene.fasta -T 45 -M 10000 -o unigene_cdhit.fasta 2</dev/null`;
print "Running TransDecoder.LongOrfs\n";
`~/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t unigene_cdhit.fasta 2</dev/null`;
$trans = "unigene_cdhit.fasta";
print "Running diamond\n";
`diamond blastp -q $trans.transdecoder_dir/longest_orfs.pep --db $uniprot_sprot_db --max-target-seqs 1 --outfmt 6 --evalue 1e-10 --threads 45 > $trans.transdecoder_dir/longest_orfs.pep.blastp.outfmt6 2</dev/null`;
print "Runing Split hmmscan\n";
`perl ~/tools/Split-Run-hmmscan.pl $trans.transdecoder_dir/longest_orfs.pep ~/PfamLib/Pfam33.1 45 $trans.transdecoder_dir/longest_orfs.pep.pfam.domtblout 2</dev/null`;
print "Runing TransDecoder.Predict\n";
`~/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t $trans --retain_pfam_hits $trans.transdecoder_dir/longest_orfs.pep.pfam.domtblout --retain_blastp_hits $trans.transdecoder_dir/longest_orfs.pep.blastp.outfmt6 2</dev/null`;
print "###########01##########\n           01  Done!!!            \nPlease See results by yourself!!!\n";


if($anno eq "yes"){
print "###########02#########\n   Run blast and hmmer   \n";
print "Runing diamond\n";
`diamond blastx --query Trinity.fasta --db $uniprot_sprot_db --threads 45 --max-target-seqs 1 --outfmt 6 > blastx.outfmt6`;
`diamond blastp --query $trans.transdecoder.pep --db $uniprot_sprot_db --threads 45 --max-target-seqs 1 --outfmt 6 > blastp.outfmt6`;

`perl ~/tools/Split-Run-hmmscan.pl Trinity.fasta.transdecoder.pep ~/PfamLib/Pfam33.1 45 TrinotatePFAM.out`;

print "###########02#########\n     02 Done!!!      \nPlease See results by yourself!!!\n";

print "###########03#########\n    Annotation Transcript   \n";
`perl ~/miniconda3/envs/Trinity/bin/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map`;
my $Trinotate_sqlite="/data/00/user/user105/Trinotate_db/Trinotate20190422.sqlite";
my $Trinotate = "~/software/Trinotate-Trinotate-v3.2.1/Trinotate";
`ln -s $Trinotate_sqlite ./Trinotate.sqlite`;
print "Runing Trinotate\nInit\n";
`$Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep`;
print "load protein hits\n";
`$Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6`;
print "load transcript hits\n";
`$Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6`;
print "load pfam hits\n";
`$Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out`;
print "report\n";
`$Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls`;
print "###########03#########\n   Annotation Finished!!!   \nCongratulations for you!!!\n";
}

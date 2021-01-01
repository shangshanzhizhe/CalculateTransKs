#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;

my ($dna,$out,$mafft) = @ARGV;
die "Usage: perl $0 <dna source file> <species name> <mafft path>\n" if(@ARGV<3);
my $rand = rand(1e5);
`mkdir -p tmp` if(! -e "tmp");
`mkdir -p tmp/temp$rand`;
my %dna_code = &DealDNAseq($dna);
`$mafft --quiet tmp/temp$rand/prot.fa > tmp/temp$rand/mafft.fa`;
my $all = Bio::SeqIO->new(-file=>"tmp/temp$rand/mafft.fa",-format=>'fasta');
open(O,">$out.best.fas") || die($!);
while(my $seq = $all->next_seq){
    my $id = $seq->id;
    my $seq = $seq->seq;
    my $length = length $seq;
    $seq =~ s/\*|\.//;
    my @a = split //,$seq;
    my @new = ();
    my $num = 0;
    foreach my $e(@a){
        if($e eq "-"){
            push(@new,"---");
        }elsif($e =~ /[A-Z]/){
            $num++;
            push(@new,$dna_code{$id}{$num});
        }
    }
    my $new_seq = join "",@new; my $len_cds = length $new_seq;
    print O ">$id\n$new_seq\n";
}
close O;
undef $all;
`rm -rf tmp/temp$rand`;
sub DealDNAseq{
    my $file = shift;
    my %hash;
    my $all = Bio::SeqIO->new(-file=>"$file",-format=>'fasta');
    open(OUT,">>tmp/temp$rand/prot.fa") || die("Can't create tmp/temp$rand/prot.fa!!!\n");
    while(my $seq = $all->next_seq){
        my $id = $seq->id;
        my $seq = $seq->seq;
        my $pep = translate_nucl($seq);
        print OUT ">$id\n$pep\n";
        my $length = length $seq;
        my $p = $length/3;
        my $mod = $length%3;
        die "The Seq length is not Multiple of 3.\n" if($mod != 0);
        for(my $i=1;$i<=$p;$i++){
            my $code = substr($seq,3*($i-1),3);
            $hash{$id}{$i} = $code;
        }
    }
    close OUT;
    return %hash;
}

sub translate_nucl{
    my $seq=shift;
    my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
    my $pro=$seq_obj->translate;
    $pro=$pro->seq;
    return($pro);
}


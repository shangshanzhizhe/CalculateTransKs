#!/usr/bin/env perl
use warnings;
use strict;

my ($gff,$out,$taxon_code) = @ARGV;
die"Usage:\nperl $0 <transdecoder gff3> <output> <taxon_code>\n" if(@ARGV < 3);

open(F,$gff) || die("Can't open the $gff!!!\n");
open(O,">$out") || die($!);
while(<F>){
    chomp;
    next if(/^#|^\s*$/);
    my @a = split /\t/,$_;
    next unless($a[2] =~ /gene/i);
    $a[8] =~ /ID=([^;]+)/;
    $1 =~ /\S+~~(\S+)/;
    print O "$a[0]\t$taxon_code\_$1\t$a[3]\t$a[4]\n";
}
close F; close O;

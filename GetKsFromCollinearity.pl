#!/usr/bin/env perl
use warnings;
use strict;

my $config = shift;
die "Usage: $0 <config file>\n" unless($config);
my %config = ReadConfig($config);
my ($out1,$lable1) = split/:/,$config{species1};
my ($out2,$lable2) = split/:/,$config{species2};
if($out1 ne $out2){
    open(O1,">$out1-$out1-Ks.out") || die($!);
    open(O2,">$out2-$out2-Ks.out") || die($!);
    open(O3,">$out1-$out2-Ks.out") || die($!);
    print O1 "Ks\tSpecies\tType\n";
    print O2 "Ks\tSpecies\tType\n";
    print O3 "Ks\tSpecies\tType\n";
    open(O,">$out1-$out2.Collinearity.Ks.out") || die($!);
    print O "Ks\tSpecies\tType\n";
    open(K,$config{collinearity}) || die("Can't open the $config{collinearity}\n");
    while(<K>){
        chomp;
        next if(/^#/);
        my @line = split/\t+/;
        #next if($line[-1] =~ /[A-Z a-z]+/);
        if($line[1] =~ /$lable1/ and $line[2] =~ /$lable1/ and $line[-1] > 0){
            print O1 "$line[-1]\t$out1-$out1\t1\n";
            print O "$line[-1]\t$out1-$out1\t1\n";
        }elsif($line[1] =~ /$lable2/ and $line[2] =~ /$lable2/ and $line[-1] > 0){
            print O2 "$line[-1]\t$out2-$out2\t1\n";
            print O "$line[-1]\t$out2-$out2\t1\n";
        }elsif($line[1] =~ /$lable1/ and $line[2] =~ /$lable2/ and $line[-1] > 0){
            print O3 "$line[-1]\t$out1-$out2\t2\n";
            print O "$line[-1]\t$out1-$out2\t2\n";
        }elsif($line[1] =~ /$lable2/ and $line[2] =~ /$lable1/ and $line[-1] > 0){
            print O3 "$line[-1]\t$out1-$out2\t2\n";
            print O "$line[-1]\t$out1-$out2\t2\n";
        }
    }
    close K; close O1; close O2; close O3; close O;
}else{
    open(O,">$out1-$out2.Collinearity.Ks.out") || die($!);
    print O "Ks\tSpecies\tType\n";
    open(K,$config{collinearity}) || die("Can't open the $config{collinearity}\n");
    while(<K>){
        chomp;
        next if(/^#/);
        my @line = split/\t+/;
        next if($line[-1] < 0);
        print O "$line[-1]\t$out1-$out2\t2\n";
    }
    close K; close O;
}

`Rscript /data/00/user/user105/tools/CalculateTransKs/Ksdistribution.R -i $out1-$out2.Collinearity.Ks.out -x 5 -o $out1-$out2.Collinearity.Ks.pdf`;

sub ReadConfig{
    my $file = shift;
    open(my $f,$file) || die("Can't open the config $file file!\n");
    my %hash;
    while(<$f>){
        chomp;
        my @a = split /\s+/;
        $hash{$a[0]} = $a[2];
    }close $f;
    return %hash;
}

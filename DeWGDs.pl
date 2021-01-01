### LZU 20200529 ###
### Created by Xin Du ###

#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
use MCE::Loop;
#use Parallel::ForkManager;
use Getopt::Long;
use POSIX;
use FindBin qw($Bin);
use Cwd qw(abs_path getcwd cwd);
use List::Util qw/sum max min/;

my $path = $Bin;
my ($input_cds,$input_pep,$config,$e_value,$thread_num,$blast,$model,$help,$pairs);
GetOptions(
           'input_cds|i=s' => \$input_cds,
           'input_pep|p=s' => \$input_pep,
           'config|c=s' => \$config,
           'model|m=s' => \$model,
           'e_value|e=s' => \$e_value,
           'thread_num|t=i' => \$thread_num,
           'blast|b=s' => \$blast,
           'cds_pairs|cp=s' => \$pairs,
           'help|?' => \$help,
          ) or die("Error in command line arguments\n");
if ( (! $config) || $help){
    &print_help;
    exit;
}
if(!$input_cds and !$pairs){
    &print_help;
    exit;
}
$model = "YN" if(! $model);
$thread_num = 4 if(! $thread_num);
$e_value = 1e-10 if(! $e_value);
my %config = ReadConfig($config);
if((!$blast) and (! $input_pep) and (!$pairs)){
    print "# No pep file or pairs dir path!!!
# The pairs dir path is prepared by manual
# The blastp results and pep file have to exists one\n";
    exit;
}
print "\nExists blastp.ok file!\nSkip blastp\n" if(-e "blastp.ok");
if(!$blast and $input_pep){
    if(! -e "blastp.ok"){
        print "\nRunning blastp!!!\n";
        RunBlast();
    }
}

my $start_time = time;
my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
print "
___Start time___
System time $now_string\n";

my (%cds,%len);
if(!$pairs){
    my ($cds,$len) = ReadSeqFile($input_cds,"cds");
    %cds = %$cds; %len = %$len;
}
`mkdir -p results` if(! -e "results" and $pairs);
`mkdir -p pairs` if($blast or $input_pep);
###prepare file###
my %class; my %inf; my %reciprocal;
my $B;
if(!$pairs){
    print "Prepare file\n";
    if($blast){
        open($B,$blast) || die("Can't open the $blast file.\n");
        %reciprocal = Pick_Up_Reciprocal_blast($blast);
    }elsif(-e "blast6.out" and !$blast){
        open($B,"blast6.out") || die("No blast6.out in the dir!!!\n");
        %reciprocal = Pick_Up_Reciprocal_blast("blast6.out");
    }else{
        die "Please you check blast results!\n";
    }
}

my %buffer; my ($sp1,$sp2)=(undef,undef);
my %pairs; my %mcl;
if(!$pairs){
    print "Deal blast result file!!!\n";
    while(<$B>){
        chomp;
        next if(/^#|^\s*$/);
        my @a = split/\s+/,$_;
        next if($a[0] eq $a[1] || $a[-2] > $e_value);
        $a[0] =~ /^(\S+)\|/; my $l1 = $1;
        $a[1] =~ /^(\S+)\|/; my $l2 = $1;
        ($l1,$l2) = sort($l1,$l2);
        die "The Seq ID is error (No | at the $a[0] or $a[1])!!!\n" unless($l1 || $l2);
        my $judge = 0;
        my $cover1 = $a[3] / $len{$a[0]};
        my $cover2 = ($a[9]-$a[8]+1) / $len{$a[1]};
        my $score = $a[2]/100 + 1 - abs($cover1-$cover2);
        $judge = 1 if($cover1 >= 0.9 and $cover2 >= 0.9);
        if(exists $reciprocal{"$a[0] $a[1]"}){
            if($l1 eq $l2 and &FilterCopy($a[0],$a[1])){
	push(@{$mcl{"$l1-$l1"}},"$a[0]\t$a[1]\t$score");
	if($a[2] >= 50 and $judge == 1){
	    $class{"$l1-$l1"}++;
	    my $class = "$l1-$l1";
	    push(@{$buffer{$class}},"$a[0]\t$a[1]");
	    my $num = sprintf("%06d",$class{$class});
	    $pairs{$class.$num}{cds} = ">$a[0]\n$cds{$a[0]}\n>$a[1]\n$cds{$a[1]}\n";
	    #$pairs{$class.$num}{pep} = ">$a[0]\n$pep{$a[0]}\n>$a[1]\n$pep{$a[1]}\n";
	}
            }elsif ($l1 ne $l2 and &FilterCopy($a[0],$a[1])){
	push(@{$mcl{"$l1-$l2"}},"$a[0]\t$a[1]\t$score");
	if($a[2] >= 50 and $judge == 1) {
	    $class{"$l1-$l2"}++;
	    ($sp1,$sp2) = ($l1,$l2);
	    my $class = "$l1-$l2";
	    push(@{$buffer{$class}},"$a[0]\t$a[1]");
	    my $num = sprintf("%06d",$class{$class});
	    $pairs{$class.$num}{cds} = ">$a[0]\n$cds{$a[0]}\n>$a[1]\n$cds{$a[1]}\n";
	    #$pairs{$class.$num}{pep} = ">$a[0]\n$pep{$a[0]}\n>$a[1]\n$pep{$a[1]}\n";
	}
            }
        }
    }
    close $B;
    die "No found that gene pairs!!!\nPlease you check input file.\n" if((keys %pairs) < 1);
    for my $type(sort keys %mcl){
        for(@{$mcl{$type}}){
            `echo "$_" >> $type.mclInuput`;
        }
    }%mcl = ();
    my $pairs_all = scalar (keys %pairs);
    my $pairs_1 = @{$buffer{"$sp1-$sp1"}};
    my $pairs_2 = @{$buffer{"$sp2-$sp2"}};
    my $pairs_3 = @{$buffer{"$sp1-$sp2"}};
    print "All found all pairs number:  $pairs_all
$sp1-$sp1 pairs number:  $pairs_1
$sp2-$sp2 pairs number:  $pairs_2
$sp1-$sp2 pairs number:  $pairs_3
Prepare have done!\n";

    %class = (); %cds = (); %len = ();
    open(O,">$sp1-$sp1.gene.pairs"); print O join"\n",@{$buffer{"$sp1-$sp1"}},"\n"; close O;
    open(O,">$sp2-$sp2.gene.pairs"); print O join"\n",@{$buffer{"$sp2-$sp2"}},"\n"; close O;
    open(O,">$sp1-$sp2.gene.pairs"); print O join"\n",@{$buffer{"$sp1-$sp2"}},"\n"; close O;
    %buffer = ();
}### prepare done ###

# parallel run
print "parallel run\n";
my @name;
if(!$pairs){
    @name = (sort keys %pairs);
}else{
    @name = <$pairs/*.cds.fa>;
    for(0..$#name){
        my @a = split/\//,$name[$_];
        $a[-1] =~ /(\S+)\.cds\.fa/;
        $name[$_] = $1;
    }
}
my $rand = sprintf("%.6f",rand(10));
"mkdir -p tmp$rand";
if($pairs){
    `echo "ID1\t\ID2\tModel\tKa\tKs\tKa\/Ks\tP_value" > results/Final_results.txt`;
}
MCE::Loop::init {chunk_size => 1,max_workers => $thread_num};
mce_loop{ &Run_Core($_) } (@name);
%pairs = ();
#read data
my %result; my %temp;
if(!$pairs){
    for (`less tmp$rand/results.txt`) {
        chomp;
        my @a = split/\t/;
        push(@{$result{$a[0]}{$a[1]}},$a[2]);
    }
    for (`less tmp$rand/temp.txt`){
        chomp;
        my @a = split/\t/;
        my $type = shift @a;
        push(@{$temp{$type}},join "\t",@a);
    }
}

###paralell done###
my $end_string = localtime; my ($d,$h,$m,$s);
$end_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
print "___Prallel End Time___\n$end_string\n";

###Print Final Results###
if(!$pairs){
    print "Print Final Results!!!\n";
    open(O1,">$sp1-$sp1.KaKs.result.txt"); open(O2,">$sp2-$sp2.KaKs.result.txt"); open(O3,">$sp1-$sp2.KaKs.result.txt");
    print O1 "ID1\tID2\tKa\tKs\tKa/Ks\tP_value\n"; print O2 "ID1\tID2\tKa\tKs\tKa/Ks\tP_value\n"; print O3 "ID1\tID2\tKa\tKs\tKa/Ks\tP_value\n";
    for my $type (keys %temp) {
        print O1 join "\n",@{$temp{$type}},"\n" if("$sp1-$sp1" eq $type);
        print O2 join "\n",@{$temp{$type}},"\n" if("$sp2-$sp2" eq $type);
        print O3 join "\n",@{$temp{$type}},"\n" if("$sp1-$sp2" eq $type);
    }
    close O1; close O2; close O3;
    %temp = ();
    print "Generate Ks.distribution.txt.\n";
    open(S,">Ks.distribution.txt"); open(ALL,">All.Type.Ks.out");
    print S "Ks\tSpecies\tType\n"; print ALL "ID\tKs_median\tKs_quartile\tKs_average\tKs_real_value\n";
    for my $name(keys %result){
        for my $type(sort keys %{$result{$name}}){
            my ($median,$quartile,$average,$real_value) = &CalculateMedian(@{$result{$name}{$type}});
            print ALL "$name\t$type\t$median\t$quartile\t$average\t$real_value\n";
            my @a=();
            @a = split/-/,$real_value if($real_value=~/-/);
            if(@a>=2){
	print S "$_\t$type\t2\n" for(@a);
            }else{
	print S "$real_value\t$type\t2\n";
            }
        }
    }close S; close ALL;
###All done###
    print "Done!!!\nTemp results file in tmp$rand dir.\nPlease see your results!!!\n";
}
my $end = localtime;
print "___Final End Time___\n",$end,"\n";
my $end_time = time;
$end_time -= $start_time; ($d,$h,$m,$s) = TransTime($end_time);
print "Sumtime is $end_time second($d day $h hour $m min $s sec)\n";
#`rm -rf pairs`;
`rm -rf tmp$rand`;
sub Run_Core{
    my $name = shift;
    if(!$pairs){
        `echo "$pairs{$name}{cds}" > pairs/$name.cds.fa`;
        `$config{prank} -d=pairs/$name.cds.fa -o=pairs/$name -once -DNA -codon -quiet 2</dev/null`;
        `perl $path/fa2axt.pl pairs/$name.best.fas pairs/$name.axt`;
        `$config{KaKs_Calculator} -i pairs/$name.axt -o pairs/$name.axt.kaks -m $model 2</dev/null`;
        &ReadKaKs("pairs/$name.axt.kaks",$name);
        `rm -f pairs/$name.cds.fa pairs/$name.axt pairs/$name.axt.kaks pairs/$name.best.fas`;
    }else{
        `$config{prank} -d=$pairs/$name.cds.fa -o=results/$name -once -DNA -codon -quiet 2</dev/null`;
        `perl $path/fa2axt.pl results/$name.best.fas results/$name.axt`;
        `$config{KaKs_Calculator} -i results/$name.axt -o results/$name.axt.kaks -m $model 2</dev/null`;
        &ReadKaKs("results/$name.axt.kaks",$name);
        `rm -f results/$name.axt results/$name.best.fas`;
    }
}

sub ReadSeqFile{
    my ($file,$type) = @_;
    my %hash;my %len;
    my $all = Bio::SeqIO->new(-file=>"$file",-format=>'fasta');
    while(my $seq = $all->next_seq){
        my $id = $seq->id;
        my $seq = $seq->seq;
        $hash{$id} = $seq;
        if($type eq "cds"){
            my $len = length($seq) / 3; $len{$id} = $len;
        }
    }
    return (\%hash,\%len) if($type eq "cds");
    return (%hash) if($type eq "pep");
}

sub ReadKaKs{
    my ($kaks,$type) = @_;
    open(F,$kaks) or die "Can't open the $kaks file!\n";
    #kaks file format
    #ID  model_name      Ka        Ks Ka/Ks        P_value
    while (<F>){
        chomp;
        next if(/^Sequence|^\s*$/);
        my @a = split/\s+/;
        next if($a[2] eq "NA" || $a[3] eq "NA" || $a[4] eq "NA" || $a[5] eq "NA" || $a[3] > 99);
        $type =~ s/\d+//g if(!$pairs);
        #print "$_\n";
        if($a[5]<0.01){
            my @id = split /;/,$a[0];
            if(!$pairs){
	`echo "$id[0]\t$type\t$a[3]" >> tmp$rand/results.txt`;
	`echo "$type\t$id[0]\t$id[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]" >> tmp$rand/temp.txt`;
            }else{
	`echo "$id[0]\t$id[1]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]" >> results/Final_results.txt`;
            }
        }
    }close F;
}

sub FilterCopy{
    my ($a1,$a2) = @_;
    if(!exists $inf{"$a1-$a2"} and !exists $inf{"$a2-$a1"}){
        $inf{"$a1-$a2"}++;
        return 1;
    }else{
        return 0;
    }
}

sub CalculateMedian{
    my @aa = @_;
    @aa = sort {$a <=> $b} @aa;
    my $first = $aa[0];
    my $hlen = @aa / 2;
    my $median = 0; my $quartile = 0; my $average = sum(@aa) / @aa;
    my $real_value;
    my (@tmp1,@tmp2,@tmp3,@tmp4,@tmp5,@tmp6);
    @tmp1 = grep {$_<=0.01} @aa;
    @tmp2 = grep {$_<=1 and $_ > 0.01} @aa;
    @tmp3 = grep {$_<=2 and $_ > 1} @aa;
    @tmp4 = grep {$_<=3 and $_ > 2} @aa;
    @tmp5 = grep {$_>3 and $_<=4} @aa;
    @tmp6 = grep {$_>4} @aa;
    if(@tmp2 or @tmp3 or @tmp4 or @tmp5 or @tmp6){
        $real_value .= "$tmp2[0]-" if($tmp2[0]);
        $real_value .= "$tmp3[0]-" if($tmp3[0]);
        $real_value .= "$tmp4[0]-" if($tmp4[0]);
        $real_value .= "$tmp5[0]-" if($tmp5[0]);
        $real_value .= "$tmp6[0]-" if($tmp6[0]);
    }else{
        $real_value = $first;
    }
    $real_value =~ s/-$//;
    if($hlen =~ /^\d+$/){
        my $tmp = $hlen/2;
        if($tmp =~ /^\d+$/){
            $quartile = ($aa[$hlen/2 -1] + $aa[$hlen/2]) / 2;
        }else{
            $quartile = $aa[$hlen/2 - 0.5];
        }
        $median = ($aa[$hlen -1] + $aa[$hlen]) / 2;
    }else{
        my $tmp = ($hlen - 0.5) / 2;
        if($tmp =~ /^\d+$/){
            $quartile = ($aa[$tmp-1] + $aa[$tmp]) / 2;
        }else{
            $quartile = $aa[$tmp - 0.5];
        }
        $median = $aa[$hlen - 0.5];
    }
    return ($median,$quartile,$average,$real_value);
}

sub RunBlast{
    `$config{diamond} makedb --in $input_pep -d database 2>/dev/null` if(!-e "blastp.ok");
    `$config{diamond} blastp --db database --query $input_pep --out blast6.out --outfmt 6 --sensitive --max-target-seqs 200 --evalue $e_value --block-size 20.0 --index-chunks 1 2</dev/null` if(!-e "blastp.ok");
    print "Blastp Done!!!\n";
    `echo "" > blastp.ok`;
}

sub ReadConfig{
    my $file = shift;
    open(my $f,$file) || die("Can't open the config $file file!\n");
    my %hash;
    while(<$f>){
        chomp;
        next if(/^#|^\s*$/);
        my @a = split /\s+/;
        $hash{$a[0]} = $a[2];
    }close $f;
    return %hash;
}

sub TransTime{
    my $sum = shift;
    my ($day,$hour,$min,$sec);
    $day = int ($sum/86400);
    $hour = int (($sum - $day*86400)/3600);
    $min = int (($sum-$day*86400 -$hour*3600)/60);
    $sec = $sum - $day*86400 - $hour*3600 - $min*60;
    return($day,$hour,$min,$sec);
}
sub Pick_Up_Reciprocal_blast{
    my $blast = shift;
    die("No all vs all blast6 result file\n") unless($blast);
    open(B,$blast) || die "Can't open the $blast!\n";
    my %inf;
    while(<B>){
        chomp;
        next if(/^#|^\s*\n$/);
        my @a = split/\s+/;
        if($a[2] >= 50){
            $inf{"$a[0] $a[1]"}++;
        }
    }

    my %reciprocal;
    for my $pair(sort keys %inf){
        my @a = split/\s+/,$pair;
        if(exists $inf{"$a[1] $a[0]"}){
            $reciprocal{"$a[0] $a[1]"}++;
        }
    }
    return %reciprocal;
}

sub print_help{
    print STDERR<<EOF;

DeWGDs (v1.00)

Usage: perl $0 --input_cds <cds_file> --input_pep <pep_file> --thread_num <thread_num> --e_value 1e-10 --config <software path config file>
       perl $0 -pairs <dir path> -t <threads number> --config|c <software path config file>
Options:
        required:
        --input_cds|i    input CDS in fasta format (The sequence name need species name lable. Example: Ath|AT2G00213)
        --input_pep|p    input PEP in fasta format (Sequence name have to same with the cds file.)
        --blast|b        option: the blastp result file(-outfmt 6)
        --config|c       software path config file
        --cds_pairs|cp   Manual generate cds pairs file in a dir, please give a path for the \"-cp\". (Notice: File name format is XXX.cds.fa)
        options:
        --model|m        Methods for estimating Ka and Ks and theirs references, Default: YN
                         NG            Nei, M. and Gojobori, T. (1986) Mol. Biol. Evol., 3, 418-426.
                         LWL           Li, W.H., Wu, C.I. and Luo, C.C. (1985) Mol. Biol. Evol., 2, 150-174.
                         LPB           Li, W.H. (1993) J. Mol. Evol., 36, 96-99.    Pamilo, P. and Bianchi, N.O. (1993) Mol. Biol. Evol., 10, 271-281.
                         MLWL          Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.
                         MLPB          Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.
                         GY            Goldman, N. and Yang, Z. (1994) Mol. Biol. Evol., 11, 725-736.
                         YN            Yang, Z. and Nielsen, R. (2000) Mol. Biol. Evol., 17, 32-43.
                         MYN           Zhang, Z., Li, J. and Yu, J. (2006) BMC Evolutionary Biology, 6, 44.
                         MS            Model Selection according to the AICc
                         MA            Model Averaging on a set of candidate models
                         GNG           Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.
                         GLWL          Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.
                         GLPB          Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.
                         GMLWL         Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.
                         GMLPB         Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.
                         GYN           Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.
                         GMYN          WangDP., Wan, HL., Zhang, S. and Yu, J. (2009) Biology Direct, 4:20 (16 June 2009)

        --e_value|e        If select --blast, don't need pick it. Default: 1e-10
        --thread_num|t     the thread number when runnig this script. Default: 4

  ** The sequence IDs within the CDS and protein files must be same!
EOF
}

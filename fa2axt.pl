use Bio::SeqIO;

my $fasta = shift;
my $out = shift;
my $all = Bio::SeqIO->new(-file=>"$fasta",-format=>'fasta');
open(O,">$out") || die($!);
my @a; my @id;
while(my $seq = $all->next_seq){
    my $id = $seq->id;
    my $seq = $seq->seq;
    push(@id,$id);
    push(@a,$seq);
}
print O $id[0],";",$id[1],"\n";
print O join "\n",@a,"\n";

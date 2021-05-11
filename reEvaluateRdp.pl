#! /usr/bin/perl -w
use strict;
use File::Which;

my $samtools = which ('samtools');

die "perl $0 <bam file>" unless (@ARGV==1 and defined $samtools);

my $poid;
my $total =0;
my $dup =0;
my $umap =0;
my $sampleID;

open IN, "$samtools view -H $ARGV[0] |" or die $!;

my @headers = <IN>;
foreach (@headers){
    if (/^\@RG/){
        $sampleID = (split /[\t\:]/)[2];
    }
    print;
}
close IN;
print "\@PG\tID:$0\tPN:$0\tCL:perl $0 $ARGV[0]\n";

open IN, "$samtools view $ARGV[0] |" or die $!;
my @item;
while (<IN>){
        $total++;
        my @a = split /\t/, $_,7;
        if ($a[1]&0x4){
            $umap++;
            print $_;
            goto ENDF;
        }
        my $po = "$a[2]\-$a[3]";
        $poid = $po if (!defined $poid);
        if ($po eq $poid ){
            push @item, \@a;
        }else{
            my $dn = &dedup (\@item);
            $poid = $po;
            $dup+=$dn;
            while (my $ele = shift @item){
                my $out = join "\t", @$ele;
                print "$out";
            }
            push @item, \@a;
        }
ENDF:
    if (eof(IN) and @item >=1){
        my $dn = &dedup (\@item);
        $dup+=$dn;
        while (my $ele = shift @item){
            my $out = join "\t", @$ele;
            print "$out";
        }
    }
}
close IN;
$sampleID = "noRG" if (!defined $sampleID);
my $dupRate = sprintf "%.4f", $dup/($total-$umap);
warn "$sampleID finds duplicate $dup in total $total with $umap ummaped\: $dupRate\n";

sub dedup {
    my $sub = shift @_;
    my %hash;
    my $dupn=0;
    foreach my $seq (@$sub){
        my $cig = $seq->[5];
        $$seq[1]-=1024 if ($$seq[1]&0x400);
        if (exists $hash{$$seq[1]} and exists $hash{$$seq[1]}{$cig}){
            $dupn++;
            $$seq[1]+=1024;
        }else{
            $hash{$$seq[1]}{$cig} =1;
        }
    }
    return ($dupn)
}

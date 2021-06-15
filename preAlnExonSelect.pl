#! /usr/bin/perl -w
use strict;

die "perl $0 <input prefix>" unless (@ARGV==1);
my $pre = shift;
my $infa = "$pre\.faa";
my $infn = "$pre\.fna";
my $inlg = "$pre\.log";
die "full list of files is not available\n" unless ( -e $infa and -e $infn and -e $inlg);
open IN, "$inlg" or die $!;
my %hash;
my %sps;
while (<IN>){
    next if (/^#/);
    chomp;
    my @a = split;
    next unless (@a==4);
    if (!exists $sps{$a[1]}){
        $sps{$a[1]} = "";
    }
    next unless ($a[3] >= 0.5);
    my @b = split /\-/, $a[2];
    next unless (@b == 2);
    my $spExon = $a[1]."-".$b[0];
    push @{$hash{$b[1]}}, $spExon;
}
close IN;
my $spn = keys %sps;
%sps = qw();
my %fna;
my %faa;
$/="\>";
open IN, "$infn" or die $!;
while (<IN>){
    chomp;
    next if ($_ eq "");
    my @a = split /\n/;
    my $sp = shift @a;
    for my $i (0..$#a){
        my $id = $sp."-".$i;
        $fna{$id} = $a[$i];
    }
}
close IN;

open IN, "$infa" or die $!;
while (<IN>){
    chomp;
    next if ($_ eq "");
    my @a = split /\n/;
    my $sp = shift @a;
    for my $i (0..$#a){
        my $id = $sp."-".$i;
        $faa{$id} = $a[$i];
    }
}
close IN;
$/="\n";
my $fnaN = keys %fna;
my $faaN = keys %faa;

die "fna ($fnaN) and faa ($faaN) contain different number of exons" unless ($fnaN == $faaN);
my $exonTaly = 0;
for my $key (sort {$a <=> $b} keys %hash){
    my @a = @{$hash{$key}};
    if ($spn == @a){
        $exonTaly++;
        foreach my $ele (@a){
            my ($sp, $num) = split /\-/, $ele;
            next unless (exists $faa{$ele});
            if (exists $sps{$sp}){
                $sps{$sp}[0].=$faa{$ele};
                $sps{$sp}[1].=$fna{$ele};
            }else{
                $sps{$sp}[0]=$faa{$ele};
                $sps{$sp}[1]=$fna{$ele};
            }            
        }
    }
}

if ($exonTaly == 0){
    die "$pre has no exon meet pre-Requirements\n";
}

my $outfa = "$pre\_ec.faa";
my $outfn = "$pre\_ec.fna";
open OFA, ">$outfa" or die $!;
open OFN, ">$outfn" or die $!;
for my $key (keys %sps){
    print OFA ">$key\n$sps{$key}[0]\n";
    print OFN ">$key\n$sps{$key}[1]\n";
}
close OFA;
close OFN;

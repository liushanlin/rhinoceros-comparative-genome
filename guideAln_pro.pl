#! /usr/bin/perl -w
use strict;
use File::Basename;

die "perl $0 <in.fasta> <sp.table>" unless (@ARGV==2);

my $faF=shift;
my $spF=shift;
my $output = basename($faF);
$output =~ s/fa/aln/ unless ($output =~ s/fasta/aln/ or $output =~ s/faa/aln/);

if (-e $output and !-z $output){
	print "$output file finished\n";
	exit (1);
}

my %spTable;
my %checkT;
open IN, "$spF" or die $!;
while(<IN>){
    my ($spQ, $spT)=(split /\s+/)[0,1];
    if ($spQ eq $spT){
        $checkT{$spQ}=1;
    }else{
		push @{$spTable{$spT}}, $spQ;
	}
}
close IN;

my $outP1=sprintf "%.4f", rand(100);
$outP1.="tmp.$output";
my %querySeqs;
open OUT, ">$outP1" or die $!;
$/="\>";
open IN, "$faF" or die $!;
while(<IN>){
    chomp;
    next if ($_ eq "");
    my ($id, $seq) = split /\n/, $_, 2;
    if (exists $checkT{$id}){
		print OUT ">$_";
        next;
    }
	$seq=~s/\n//g;
	$querySeqs{$id}=$seq;    
}
close IN;
close OUT;
my $outP2="$outP1\.aln";
my $judge = system ("mafft --anysymbol --auto --quiet $outP1 > $outP2");

die "mistakes happend during mafft alignment\n" unless (defined $judge);

open IN, "$outP2" or die $!;
open OUT, ">$output" or die $!;
while (<IN>){
    chomp;
    next if ($_ eq "");
    my ($id, $seq) = split /\n/, $_, 2;
	$seq =~s/\n//g; 
    unless (exists $spTable{$id}){
        print OUT ">$id\n$seq\n";
        next;
    }
    $seq = uc($seq);
    my $baseNum = $seq;
	$baseNum =~ s/-//g;
	$baseNum = length $baseNum;
    my @seqAray = split /\s*/, $seq;
    TMD:foreach my $ele (@{$spTable{$id}}){
		next TMD unless (exists  $querySeqs{$ele});
        my $querySeq = $querySeqs{$ele};
        my $qSL = length ($querySeq);
        unless ($qSL == $baseNum){
            unlink("$outP1") or die "cannot remove $outP1";
            unlink("$outP2") or die "cannot remove $outP2";
            die "$faF\tseq length is differ for $id and $ele";
        }
        my @querySeqArray = split /\s*/, $querySeq;
        my $j=0;
        for my $i (0..$#seqAray){
            if ($seqAray[$i] ne "-"){
                $seqAray[$i] = $querySeqArray[$j];
                $j++;
            }
        }
        my $alnedSeq = join "", @seqAray;
        print OUT ">$ele\n$alnedSeq\n";
    }
    print OUT ">$id\n$seq\n";
}
close IN;
close OUT;
unlink("$outP1") or die "cannot remove $outP1";
unlink("$outP2") or die "cannot remove $outP2";

#! /usr/bin/perl -w
use strict;
use File::Basename;

die "perl $0  <exonerate.log> <genes Maf file> <out prefix>" unless (@ARGV==3);

my $f1=shift;
my $f2=shift;
my $out=shift;

my %seqSet;
open IN, "$f1" or die $!;
while (<IN>){
    chomp;
    next unless (/vulgar/);
    my @a=split /\s+/, $_;
    my ($qp,$qs,$qe,$sp,$sps,$spe,$spd)=@a[1,2,3,5,6,7,8];
    my @ordert;
    my @orderq;
    my @tags;
    my $total;
    for (my $i=10; $i<=$#a-2; $i+=3){
        my ($label, $query, $target)=@a[$i, $i+1, $i+2];
        die "$_ is not well splited" unless ($label eq "M" or $label eq "G" or $label eq "5" or $label eq "3" or $label eq "I" or $label eq "S" or $label eq "F");
        push @ordert, $target;
        push @orderq, $query;
        push @tags, $label;
        $total+=$target;
    }
    my $length = abs($sps - $spe);
    die "$length\t$total\t$sp\t cannot get all the elements\n" unless ($length == $total);
    my @framePosition = &getFramePos(\$qs, \@orderq, \$sp, \@ordert, \$sps, \$spe, \$spd, \@tags);
    my $framePositions = join "-", @framePosition;
    if (exists $seqSet{$sp}){
        warn "$f1 file contains duplication\n";
        exit(1);
    }else{
        $seqSet{$sp}= $framePositions;
    }
}
close IN;

open IN, "$f2" or die $!;
my %start;
my %spPos;
$/="";
while (my $mafFile=<IN>){
   chomp($mafFile);
   my @lines = split /\n/, $mafFile;
   TTT:foreach my $ele (@lines){
        next TTT unless ($ele=~/^s/);
        my @a=split /\s+/, $ele,7;
        my ($spm, $chr)=(split /\./,$a[1],2)[0,1];
        my ($st,$len,$di,$mafLen)=@a[2,3,4,5];
        if($di eq "-"){
            $st=$mafLen-$st-$len;
        }
        if (exists $start{$spm}){
            $start{$spm}+=$len;
        }else{
            $start{$spm}=0;
        }
        next TTT unless (exists $seqSet{$spm});
        my @frames = split /\-/, $seqSet{$spm};
        TNT: for (my $i =0; $i <= $#frames; $i++){
            next TNT if (defined $spPos{$spm}[$i]);
            if ($frames[$i] <= $start{$spm}){
                my $tradeoff = $start{$spm}-$frames[$i]+$st;
                print "$spm\t$start{$spm}\t$frames[$i]\t$st\t$tradeoff\n";
                my $info = $chr.":".$tradeoff;
                $spPos{$spm}[$i] = $info;
            }
        }
    }
}
close IN;
$/="\n";

open LGF, ">$out\.log" or die $!;
for my $key (keys %spPos){
    for my $ele (@{$spPos{$key}}){
        print LGF "$key\t$ele\n"
    }
}
close LGF;


sub getFramePos {
    my ($p1, $p2, $p3, $p4, $p5, $p6, $p7, $p8)=@_[0,1,2,3,4,5,6,7];
    my $qs=$$p1;
    my @orderq = @$p2;
    my $sp = $$p3;
    my @ordert = @$p4;
    my $sps = $$p5;
    my $spe = $$p6;
    my $spd = $$p7;
    my @tags = @$p8;
    my @fpos;
    if ($spd eq "+"){
        for my $i (0..$#ordert){
            my $ele = $ordert[$i];
            my $tag = $tags[$i];
            my $ql = $orderq[$i];
            $sps+=$ele;
            if ($tag eq "F"){
                push @fpos, $sps;
            }
        }
    }elsif ($spd eq "-"){
       $sps=$spe;
       @ordert = reverse @ordert;
       @tags = reverse @tags;
       @orderq = reverse @orderq;
       for my $i (0..$#ordert){
            my $ele = $ordert[$i];
            my $tag = $tags[$i];
            my $ql = $orderq[$i];
            $sps+=$ele;
            if ($tag eq "F"){
                push @fpos, $sps;
            }
        }
    }else{
        die "$spd is not formated as + or -\n"
    }
    return @fpos;
}

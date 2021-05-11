#! /usr/bin/perl -w
use strict;
use lib "/home/liushanlin/testRhino/command";
use MyAlignSM qw(poscore score);
use Bio::Seq qw(new;translate;seq);
use File::Basename;

die "perl $0  <query faa |exon per line> <sp table> <exonerate.log> <input target fasta> <input query fasta> <out prefix>" unless (@ARGV==6);

my $queryfaa =shift;
my $sptable =shift;
my $f1=shift;
my $f2=shift;
my $f3=shift;
my $out=shift;



my @qfas;
my $qfa;
my @qlen;
my $qtl =0;
open IN, "$queryfaa" or die $!;
while(<IN>){
    chomp;
    next if (/^>/ or $_ eq "");
    my $len = length $_;
    push @qfas, $_;
    $qfa.=$_;
    push @qlen, $len;
    $qtl += $len;
}
close IN;
my %four;
for (my $i=0; $i <= $qtl-4; $i++){
    my $fourmer = substr $qfa, $i, 4;
    push @{$four{$fourmer}}, $i;
}

my %spInfo;
open IN, "$sptable" or die $!;
while(<IN>){
    chomp;
    my @a=split;
    die "cannot splited into 3 for $_\n" unless (@a==3);
    next if ($a[0] eq $a[1]);
    push @{$spInfo{$a[1]}},$a[0];
}
close IN;

$/="\>";
my %seqs;
open IN, "$f2" or die $!;
while (<IN>){
    chomp;
    next if ($_ eq "");
    my ($id, $seq)=(split /\n/, $_, 2)[0,1];
    $seq=~s/\n//g;
    $seqs{$id}=$seq;
}
close IN;

open IN, "$f3" or die $!;
while (<IN>){
    chomp;
    next if ($_ eq "");
    my ($id, $seq)=(split /\n/, $_, 2)[0,1];
    $seq=~s/\n//g;
    $seqs{$id}=$seq;
}
close IN;
$/="\n";

my %seqSet;
my %exonPos;
open IN, "$f1" or die $!;
while (<IN>){
    chomp;
    next unless (/vulgar/);
    my @a=split /\s+/, $_;
    my ($qp,$qs,$qe,$sp,$sps,$spe,$spd)=@a[1,2,3,5,6,7,8];
    next, warn "$sp doesn't exists in the fasta file" unless (exists $seqs{$sp});
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
    my $seqGe = &getExon(\$qs, \@orderq, \$sp, \@ordert, \$sps, \$spe, \$spd, \@tags);
    if (exists $seqSet{$sp}){
        warn "$f1 file contains duplication\n";
        exit(1);
    }else{
        $seqSet{$sp}= $seqGe;
    }
    if (exists $spInfo{$sp}){
        foreach my $sptarget (@{$spInfo{$sp}}){
            my $seqGe = &getExon(\$qs, \@orderq, \$sptarget, \@ordert, \$sps, \$spe, \$spd, \@tags);
            $seqSet{$sptarget}= $seqGe;
        }
    }
}
close IN;

open DNA, ">$out\.fna" or die $!;
open PRO, ">$out\.faa" or die $!;
open LGF, ">$out\.log" or die $!;
my $exonNum;
for my $key (keys %seqSet){
    my @seqs = split /\+/, $seqSet{$key};
    $exonNum = $#seqs if (!defined $exonNum);
    warn "in the $f1 file $key has exon of $#seqs, not equl to $exonNum\n" if ($exonNum != $#seqs);
    print DNA ">$key\n";
    print PRO ">$key\n";
    foreach my $i (0..$#seqs){
        if ($seqs[$i]=~s/\#//g){
            print LGF "#framshift\t$key\t$i\n";
        }
        my $seq=Bio::Seq -> new (-seq => $seqs[$i]);
        my $pro=$seq -> translate();
        my $faa=$pro->seq();
        my $alnS = &getPos ($faa);
        my $dis;
        if ($alnS == 1000000){
            ($alnS, $dis) = &poscore($qfa, $faa);
        }else{
            my $subread = substr $qfa, $alnS, length($faa);
            $dis = &score($subread, $faa);
        }
        $dis = sprintf "%.2f", $dis;
        my $refExon = &exonFind($alnS);
        print LGF "exonScore\t$key\t$i\-$refExon\t$dis\n";
        print DNA "$seqs[$i]\n";
        print PRO "$faa\n";
    }
}
close DNA;
close PRO;
close LGF;

sub getPos {
    my $seq = shift @_;
    my $sl = length $seq;
    my %sts;
    PTT:for (my $i =0; $i <= $sl-4; $i++){
        my $fmer = substr ($seq, $i, 4);
        next unless (exists $four{$fmer});
        foreach my $ele (@{$four{$fmer}}){
            my $st = $ele - $i;
            if (exists $sts{$st}){
                $sts{$st}++;
            }else{
                $sts{$st}=1;
            }
        }
    }
    my @out;
    my $c=0;
    PTN:for my $key (sort {$sts{$b} <=> $sts{$a}} keys %sts){
        push @out, $key;
        $c++;
        last PTN if ($c>=2);
    }
    if (@out == 1){
        return $out[0];
    }elsif (@out >=2 and $sts{$out[0]} - $sts{$out[1]} >= 2){
        return $out[0];
    }else{
        return 1000000;
    }
}


sub getExon {
    my ($p1, $p2, $p3, $p4, $p5, $p6, $p7, $p8)=@_[0,1,2,3,4,5,6,7];
    my $qs=$$p1;
    my @orderq = @$p2;
    my $sp = $$p3;
    my @ordert = @$p4;
    my $sps = $$p5;
    my $spe = $$p6;
    my $spd = $$p7;
    my @tags = @$p8;
    my $seq = uc($seqs{$sp});
    my $spSeq="";
    my $scount =0;
    my $icount =0;
    if ($spd eq "+"){
        for my $i (0..$#ordert){
            my $ele = $ordert[$i];
            my $tag = $tags[$i];
            my $ql = $orderq[$i];
            if ($ql == 0 ){
                if ($tag eq "S"){
                    $scount++;
                    $spSeq .= substr $seq, $sps, $ele;
                }
                if ($tag eq "I" and $scount == 0){
                    $icount++;
                    $spSeq.="+";
                }
                if ($tag eq "F"){
                    for my $n (1..$ele){
                        $spSeq.="#";
                    }
                }
                $sps+=$ele;
            }elsif($ql > 0){
                $spSeq .= substr $seq, $sps, $ele;
                $sps+=$ele;
                if ($tag eq "S"){
                    $scount--;
                    $spSeq.="+";
                }
            }else{
                die "$ql value does not fit\n"
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
            if ($ql == 0 ){
                if ($tag eq "S"){
                    $scount--;
                    $spSeq .= substr $seq, $sps, $ele;
                }
                if ($tag eq "I" and $scount == 0){
                    $icount++;
                    $spSeq.="+";
                }
                if ($tag eq "F"){
                    for my $n (1..$ele){
                        $spSeq.="#";
                    }
                }
                $sps+=$ele;
            }elsif($ql > 0){
                if ($tag eq "S"){
                    $scount++;
                    $spSeq.="+";
                }
                $spSeq .= substr $seq, $sps, $ele;
                $sps+=$ele;
            }else{
                die "$ql value does not fit\n"
            }
        }
        $spSeq =~ tr/ATCG/TAGC/;
        $spSeq = reverse $spSeq;
    }else{
        die "$spd is not formated as + or -\n"
    }
    return $spSeq;
}


sub getSeq {
    my ($p1, $p2, $p3, $p4, $p5)=@_[0,1,2,3,4];
    my $sp=$$p1;
    my @order = @$p2;
    my $sps = $$p3;
    my $spe = $$p4;
    my $spd = $$p5;
    my $seq = uc($seqs{$sp});
    my $spSeq="";
    if ($spd eq "+"){
        foreach my $ele (@order){
            if ($ele=~s/\+$//){
                $spSeq .= substr $seq, $sps, $ele;
                $sps+=$ele;
            }elsif ($ele =~ s/\-$//){
                $sps+=$ele;
            }else{
                die "$ele not formated with + or - at the end \n"
            }
        }
    }elsif ($spd eq "-"){
        $sps=$spe;
        foreach my $ele (reverse @order){
            if ($ele=~s/\+$//){
                $spSeq .= substr $seq, $sps, $ele;
                $sps+=$ele;
            }elsif ($ele =~ s/\-$//){
                $sps+=$ele;
            }else{
                die "$ele not formated with + or - at the end \n"
            }
        }
        $spSeq =~ tr/ATCG/TAGC/;
        $spSeq = reverse $spSeq;
    }else{
        die "$spd is not formated as + or -\n"
    }
    return $spSeq;
}

sub exonFind {
    my $po = shift @_;
    my $lenAcu = 0;
    my $e = 0;
    TNT:for my $i (0..$#qlen){
        $lenAcu+=$qlen[$i];
        if ($po <= $lenAcu){
            my $r = $lenAcu - $po;
            my $ratio = $r/$qlen[$i];
            $e = $ratio >= 0.5 ? $i : $i+1;
            last TNT;
        }
    }
    return $e;
}

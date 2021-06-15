#! /usr/bin/perl -w
use strict;
use Bio::Seq qw(new;translate;seq);

die "perl $0 <alignment> <output> <strict merge or not, 1 or 0>" unless (@ARGV==3);

my $in = shift;
my $out = shift;
my $strict = shift;

open IN, "$in" or die $!;

my %hash;
my $len;
$/="\>";
while(<IN>){
    chomp;
    next if ($_ eq "");
    my ($id, $seq) = (split /\n/, $_, 2)[0,1];
    $seq=~s/\n//g;
    my $leng = length $seq;
    $len = $leng if (!defined $len);
    die "for file $in, $id has a length of $leng, differ from the expected length of $len\n" unless ($len == $leng);
    my $uid = substr ($id, 0, 5);
    push @{$hash{$uid}}, $seq;
}
close IN;

open OUT, ">$out" or die $!;
for my $key (keys %hash){
    my @a = @{$hash{$key}};
    if (@a ==1){
        print OUT ">$key\n$a[0]\n";
    }elsif (@a == 2){
        my $seq;
        if ($strict == 1){
            $seq = &mergeReadStrict(\@a);
        }else{
            $seq = &mergeRead (\@a);
        }
        warn "$key does not merge sucessfully" unless (defined $seq);
        print OUT ">$key\n$seq\n";
    }
}

sub mergeRead {
    my $arrary = shift @_;
    my $seqf = shift @$arrary;
    my $seqs = shift @$arrary;

    my $fnaf=Bio::Seq -> new (-seq => $seqf);
    my $fnas=Bio::Seq -> new (-seq => $seqs);
    my $prof=$fnaf -> translate(-terminator => 'X');
    my $pros=$fnas -> translate(-terminator => 'X');
    my $faaf=$prof->seq();
    my $faas=$pros->seq();


    my @pf = split /\s*/, $faaf;
    my @ps = split /\s*/, $faas;
    my @df = split /\s*/, $seqf;
    my @ds = split /\s*/, $seqs;

    my $finalSeq;
    for my $i (0..$#pf){
        my $j=$i*3;
        my $subf = substr($seqf,$j,3);
        my $subs = substr($seqs,$j,3);
        if ($subf eq $subs){
            $finalSeq.=$subf;
        }else{
            my $ch;
            if ($pf[$i] ne "X" and $ps[$i] ne "X"){
                $ch = "NNN"
            }elsif ($pf[$i] eq "X"){
                $ch = $subs;
            }else{
                $ch = $subf;
            }
            $finalSeq.=$ch;
        }
    }
    return $finalSeq;
}

sub mergeReadStrict {
    my $arrary = shift @_;
    my $seqf = shift @$arrary;
    my $seqs = shift @$arrary;

    my @df = split /\s*/, $seqf;
    my @ds = split /\s*/, $seqs;

    my $finalSeq;
    for(my $i=0; $i<$#df;$i+=3){
        my $subf = substr($seqf,$i,3);
        my $subs = substr($seqs,$i,3);
        if ($subf eq $subs){
            $finalSeq.=$subf;
        }else{
            $finalSeq.="NNN";
        }
    }
    return $finalSeq;
}

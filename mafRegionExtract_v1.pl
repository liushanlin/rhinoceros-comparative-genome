#! /usr/bin/perl -w
use strict;
use File::Basename;


die "perl $0 <maf input> <target species list> <maximum gap to be concatinated | 50 by default>" unless (@ARGV==3);

my $maf=shift;
my $spt=shift;
my $cut=shift;
$cut = defined $cut ? $cut:50;
my $spnum=0;
my %sps;
my @species;
open IN, "$spt" or die $!;
while(<IN>){
    chomp;
    my $spn=(split /\s+/)[0];
    $sps{$spn}=1;
    $spnum++;
    push @species, $spn;
}
close IN;
my %info;
my $block=0;
open IN, "$maf" or die $!;
$/="";
TNT:while (<IN>){
    chomp;
    my $record;
    my $judge=0;
    my @lines = split /\n/;
    TTT:foreach my $ele (@lines){
        next TTT unless ($ele=~/^s/);
        my @a=split /\s+/, $ele,7;
        my ($spm, $chr)=(split /\./,$a[1],2)[0,1];
        next unless (exists $sps{$spm} and defined $chr);
        $judge++;
        my ($st,$len,$di,$mafLen)=@a[2,3,4,5];
        my $end;
        if ($di eq "+"){
            $end=$st+$len;
            $st++;
        }elsif($di eq "-"){
            $end=$mafLen-$st;
            $st=$mafLen-$st-$len;
        }else{
            die "$di is not + or - \n";
        }       
        my $pos="$spm\t$chr:$st\:$end\:$di\n";
        $record.=$pos;
    }
    next TNT unless ($judge == $spnum);
    foreach my $ele (split /\n/, $record){
        my ($spm, $pos)=(split /\t/, $ele)[0,1];
        push @{$info{$spm}},$pos;
    }
    $block++;
}
close IN;
$/="\n";

my %candidate;
TNN:for my $i (0..$block-2){
    my $j=$i+1;
    foreach my $sp (@species){
        my @a=split /\:/, $info{$sp}[$i];
        my @b=split /\:/, $info{$sp}[$j];
        my ($f_chr,$f_s,$f_e,$f_d)=(split /\:/, $info{$sp}[$i])[0,1,2,3];
        my ($r_chr,$r_s,$r_e,$r_d)=(split /\:/, $info{$sp}[$j])[0,1,2,3];
        if ($f_chr eq $r_chr and $f_d eq $r_d){
            if ($f_d eq "+"){
				next TNN unless (abs($r_s-$f_e)<=$cut);
			}else{
				next TNN unless (abs($f_s-$r_e)<=$cut);
			}
        }else{
            next TNN;
        }
    }
    $candidate{$j}=$i;
}
my @region;
for my $i (0..$block-1){
    push @region, $i unless (exists $candidate{$i});
}

for my $i (0..$#region-1){
    my $j=$i+1;
    my $out;
    if ($region[$j]-$region[$i]==1){
        $out = &printBlock($region[$i]);
        print ">$region[$i]\n$out";
    }else{
        my $last=$region[$j]-1;
        my $putin = "$region[$i]\-$last";
        $out= &printBlock($putin);
        print ">$putin\n$out";
    }
    if($j==$#region){
        $out = &printBlock($region[$j]);
        print ">$region[$j]\n$out";
    }
}

exit(1);

sub printBlock {
    my $bn=shift @_;
    my @a=split /\-/, $bn;
    my $outarg;
    foreach my $sp (@species){
        if (@a==1){
            my ($f_chr,$f_s,$f_e,$f_d)=(split /\:/, $info{$sp}[$bn])[0,1,2,3];
            $outarg.="$sp\t$f_chr\:$f_s\-$f_e\t$f_d\n";
        }elsif(@a==2){
            my ($a1, $a2)=@a[0,1];
            my ($f_chr,$f_s,$f_e,$f_d)=(split /\:/, $info{$sp}[$a1])[0,1,2,3];
            my ($r_chr,$r_s,$r_e,$r_d)=(split /\:/, $info{$sp}[$a2])[0,1,2,3];
            if ($f_d eq $r_d and $f_d eq "+"){
                $outarg.="$sp\t$f_chr\:$f_s\-$r_e\t$f_d\n"; 
            }elsif($f_d eq $r_d and $f_d eq "-"){
                $outarg.="$sp\t$f_chr\:$r_s\-$f_e\t$f_d\n";
            }else{
                die "$info{$sp}[$a1] and $info{$sp}[$a2] direction conflict";
            }
        }else{
            die "$bn cannot be separated to 2 elements by -\n";
        }
    }
    return $outarg;
}

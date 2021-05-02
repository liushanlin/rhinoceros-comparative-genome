#! /usr/bin/perl -w
use strict;

die "perl $0 <sp.table> <ref chr.info> <blocks> <window size | 100K by default>" unless (@ARGV>=3);

my $spt=shift;
my $chrf=shift;
my $bf=shift;
my $size = @ARGV==4? shift:100000;

my %reference;
my %spInfo;
my @spOrder;

open IN, "$spt" or die $!;
while(<IN>){
    chomp;
    my @a=split;
    die "cannot splited into 3 for $_\n" unless (@a==3);
    die "dup $a[0]\n" if (exists $spInfo{$a[0]});
    @{$spInfo{$a[0]}}[0,1]=@a[1,2];
    push @spOrder, $a[0];
}
close IN;

my %chrs;
open IN, "$chrf" or die $!;
while(<IN>){
    chomp;
    my @a=split;
    die "cannot splited into 2 for $_\n" unless (@a==2);
    die "dup $a[0]\n" if (exists $chrs{$a[0]});
    $chrs{$a[0]}=$a[1];
}
close IN;

open IN, "$bf" or die $!;
<IN>;
my $refF=<IN>;
my $ref=(split /\s+/, $refF)[0];
warn "$ref is the reference species, isn't it? 
if not, please stop the prgram\nif yes, ignore this message\n";

sleep(10);

my %coord;
seek(IN,0,0);
$/="\>";
<IN>;
my $spNum;
TNT:while(<IN>){
    chomp;
    my @lines = split /\n/, $_, 3;
    $spNum = $#lines if (!defined $spNum);
    my $blockNum = shift @lines;
    warn "sp num in $blockNum is $#lines, not equall to $spNum\n" unless ($spNum == $#lines+1);
    my $ele = shift @lines;
    my ($sp, $chrPos, $ori) = (split /\s+/, $ele)[0,1,2];
    die "$ele is not correctly splited" unless (defined $ori and ($ori eq "+" or $ori eq "-") and $sp eq $ref);
    my ($chr,$start)=(split /[\:\-]/, $chrPos)[0,1];
    die "$chrPos is not correclty splited\n" unless (defined $chr and defined $start);
    next TNT unless (exists $chrs{$chr});
    my $hChr=$chrs{$chr};
    my $sizeOrder = int($start/$size)+1;
    push @{$coord{$hChr}{$sizeOrder}},$_;
}
close IN;

open TMP, ">tmpFile1" or die $!;
for my $key (keys %coord){
    for my $sk (sort {$a <=> $b} keys %{$coord{$key}}){
        my $region=$sk*$size;
        my $out="$key\_$region";
        print TMP ">$out\n";
        my @blockInfo = @spOrder;
        foreach my $ele (@{$coord{$key}{$sk}}){
            my @spPos = &getInfo($ele);
            warn "number conflict for $key $sk $ele blocks" unless ($#spPos == $#blockInfo);
            for my $i (0..$#spPos){
                $blockInfo[$i].="\t$spPos[$i]";
            }
        }
        foreach my $ele (@blockInfo){
            print TMP "$ele\n";
        }
    }
}
close TMP;

%coord=();

print "coordination file was generated and going to extract sequences next
you can stop here by type Ctr+C if you don't want to sequences\n";

sleep (10);

my %chrCheck;
for my $key (keys %reference){
    my $id = (split /\:/, $key, 2)[0];
    push @{$chrCheck{$id}},$key;
}

%reference = ();

my %bNumC;

foreach my $spEle (@spOrder){
    my $genome = $spInfo{$spEle}[1];
    my %genoSeqs;
    open FA, "$genome" or die $!;
    TMD:while (my $line=<FA>){
        chomp($line);
        next TMD if ($line eq "");
        my ($sid, $seq)=(split /\n/, $line,2)[0,1];
        $seq=~s/\n//g;
        my $chrName = (split /\s+/, $sid, 2)[0];
        my $check = "$spEle\_$chrName";
        next TMD unless (exists $chrCheck{$check});
        foreach my $block (@{$chrCheck{$check}}){
            my ($posSE, $ori) = (split /\:/, $block)[-2,-1];
            my ($start, $end) = (split /\-/, $posSE)[0,1];
            die "not well formated $block\ngenerate $start\t$end\t$ori\n" unless (defined $start and defined $end and ($ori eq "+" or $ori eq "-"));
            my $len=$end-$start+1;
            my $subSeq = substr $seq, $start, $len;
            if ($ori eq "-"){
                $subSeq = uc($subSeq);
                $subSeq=~tr/ATCGN/TAGCN/;
                $subSeq = reverse ($subSeq);
            }
            $genoSeqs{$block}=$subSeq;
        }
    }
    close FA;
    
    warn "$spEle has been loaded sucessfully\n";

    open IN, "tmpFile1" or die $!;
    TME:while(my $bl=<IN>){
        chomp($bl);
        next TME if ($bl eq "");
        my @lines = split /\n/, $bl;
        my $id = shift @lines;
        my $chrNum = (split /\_/, $id)[0];
        mkdir $chrNum unless ( -d $chrNum);
        open OUT, ">>$chrNum/$id\.fasta" or die $!;
        my $eleNum;
        TMF:foreach my $ele (@lines){
            my @a = split /\t/, $ele;
            if (exists $bNumC{$id}){
                die "$id blocks doesn't match\n" unless ($#a == $bNumC{$id})
            }else{
                $bNumC{$id}=$#a;
            }
            my $spName = shift @a;
            next TMF unless ($spName eq $spEle);
            my @tmpSeqs;
            foreach my $blockID (@a){
                my $keyC="$spName\_$blockID";
                die "$keyC don't have data\n" unless (exists $genoSeqs{$keyC});
                push @tmpSeqs, $genoSeqs{$keyC};
            }
            my $mergedSeq = join "", @tmpSeqs;
            print OUT ">$spName\n$mergedSeq\n";
        }
        close OUT;
    }
    close IN;
    %genoSeqs=();
}

%chrCheck=();

exit(0);

sub getInfo {
    my $block = shift @_;
    my @lines = split /\n/, $block;
    shift @lines;
    my %sinfo;
    foreach my $ele (@lines){
        my ($sp, $chrPos, $ori) = (split /\s+/, $ele)[0,1,2];
        $sinfo{$sp}="$chrPos\:$ori";
    }
    my @seqs;
    foreach my $ele (@spOrder){
        my $target = $spInfo{$ele}[0];
        my $id = "$ele\_$sinfo{$target}";
        warn "$id dup when record reference hash\n" if (exists $reference{$id});
        $reference{$id}=1;
        push @seqs, $sinfo{$target};
    }
    return @seqs;
}

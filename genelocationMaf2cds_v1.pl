#! /usr/bin/perl -w
use strict;
use File::Basename;
use Bio::Seq qw(new translate seq);
use Getopt::Long;

=head1 usage

	perl $0 <parameters>
	
		-gen 	genelocaiton file
		-maf 	list of maf files
		-ref 	ref species in the maf file
		-cut 	cutoff of gene integrity | 0.8 by default
		-min 	minimum length of Amino acid into the output file | 10 bp by default
		-ski    will skip CDSs which have stop codon in between (useful for bad annotated genome reference)
		-help 	print out current help info

=cut


GetOptions(
        "gen:s"=>\our$f,
        "maf:s"=>\our$maf,
        "ref:s"=>\our$spt,
        "cut=i"=>\our$cut,
        "min=i"=>\our$minP,
		"ski"=>\our$Ski,
        "help"=>\our$Help
);

die `pod2text $0` if (!defined $f or !defined $maf or !defined $spt or defined $Help);

$cut = defined $cut ? $cut:0.8;
$minP = defined $minP ? $minP:10;

warn "CDSs which have stop codon in between will be skipped\n" if (defined $Ski);

my %gl;

open IN, "$f" or die $!;
$/="gene";
while (<IN>){
	chomp;
	next if ($_ eq "");
	my @a=split /\n/;
	my $id = (split /\t/,$a[0])[0];
	my %hash;
	my $strand;
	TNT:for my $i (1..$#a){
		next TNT if ($a[$i] eq "");
		my ($st,$en,$fr,$di)=(split /\t/,$a[$i])[0,1,2,3];
		$strand = $di if (!defined $strand);
		$fr = 0 if ($fr eq '.');
		warn "not same direction for gene $id at start pos $st" unless ($strand eq $di);
		if (exists $hash{$st}){
			if ($en > $hash{$st}){
				$hash{$st}[0]=$en;
				$hash{$st}[1]=$fr;
			}
		}else{
			$hash{$st}[0]=$en;
			$hash{$st}[1]=$fr;
		}
	}
	push @{$gl{$id}},$strand;
	for my $key (sort keys %hash){
		push @{$gl{$id}},($key,$hash{$key}[0],$hash{$key}[1]);
	}
	undef %hash;
}
close IN;
$/="\n";
open IN, "$maf" or die $!;
TOL:while (<IN>){
	chomp;
	my $mf=basename($_);
	my $id=(split /\./,$mf)[0];
	$id=~s/^gene//;
	my %seq;
	my %cds;
	my $tl;
	unless (exists $gl{$id}){
		warn "why no $id in the genelocation file\n";
		next TOL;
	}
	for (my $i=1;$i<$#{$gl{$id}};$i+=3 ){
		my $n=$i+1;
		my $nn=$i+2;
		my ($st,$en,$fr)=($gl{$id}[$i],$gl{$id}[$n],$gl{$id}[$nn]);
		$seq{$st}[0]=$en;
		$seq{$st}[1]=$fr;
		my $cl=$en-$st+1;
		$tl+=$cl;
	}
	my $gdi=$gl{$id}[0];
	my $last=0;
	my @part;
	open MF, "$_" or die $!;
	TTT:while(<MF>){
		next TTT unless (/^s/);
		my @a=split;
		next TTT unless ($a[1]=~/^$spt/);
		my ($st,$len,$di,$seq)=@a[2,3,4,6];
		if ($st==$last){
			$part[-1].=$seq;
		}else{
			push @part, ($st,$di,$seq);
		}
		$last=$st+$len;
	}
	close MF;
	next TOL if ($last==0);
	my $cdsLen=0;
	for (my $i=0;$i<$#part;$i+=3){
		my $ii=$i+1;
		my $iii=$i+2;
		my ($st,$di,$seq)=@part[$i,$ii,$iii];
		$seq=~s/\-//g;
		my $len=length $seq;
		my $en=$st+$len;
		for my $key (keys %seq){
			if ($st<=$key and $en>=$seq{$key}[0]){
				my $cst=$key-$st-1;
				my $clen=$seq{$key}[0]-$key+1;
				while ($cst < 0 ) {
					$cst+=3;
					$clen-=3;
				}
				$cdsLen+=$clen;
				my $s=substr ($seq,$cst,$clen);
				if ($di eq $gdi){
					$cds{$key}[0]=$s;
					$cds{$key}[1]=$seq{$key}[1];
				}else{
					$s=uc $s;
					$s=~tr/ATCGN/TAGCN/;
					$s=reverse $s;
					$cds{$key}[0]=$s;
					$cds{$key}[1]=$seq{$key}[1];
				}
			}
		}
	}
	my $cdsNum=keys %cds;
	next TOL unless ($cdsNum >= 1);
	my $extR=$cdsLen/$tl;
	next TOL unless ($extR >= $cut);
	open PRO, ">gene$id\.faa" or die $!;
	print PRO ">gene$id\t$gdi\n";
	if ($gdi eq "+"){
		TNN:for my $key (sort keys %cds){
			next TNN if ((length $cds{$key}[0])< 3);
			my $fna=Bio::Seq -> new (-seq => $cds{$key}[0]);
			my $faa= $fna -> translate(-frame => $cds{$key}[1]);
			my $pro= $faa -> seq();
			my $tmpP=$pro;
			$tmpP=~s/\*$//;
			my $stop = $tmpP=~s/\*//g;
			next TNN if ($stop >= 1 and defined $Ski);
			print PRO "$pro\n" if ((length $pro) >= $minP);
		}
	}else{
		NNN:for my $key (sort {$b <=> $a}keys %cds){
			next NNN if ((length $cds{$key}[0])< 3);
			my $fna=Bio::Seq -> new (-seq => $cds{$key}[0]);
			my $faa= $fna -> translate(-frame => $cds{$key}[1]);
			my $pro= $faa -> seq();
			my $tmpP=$pro;
			$tmpP=~s/\*$//;
			my $stop = $tmpP=~s/\*//g;
			next NNN if ($stop >= 1 and defined $Ski);
			print PRO "$pro\n" if ((length $pro) >= $minP);			
		}		
	}
	close PRO;
}

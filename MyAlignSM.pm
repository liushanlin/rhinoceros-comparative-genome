#! MyAlignSM.pm
use strict;
package MyAlignSM;
use Exporter qw(import);
our @EXPORT_OK = qw (single twice poscore score);

my $Mat = 2;
my $Mis = 1;
my $Gap = 5;
my $Gex = 1;
#my %jud={"A_A",1,"T_T",1,"C_C",1,"G_G",1}; #define your own score
my (@Sco, @Poi, @Poj);
my $t_max=0;
my ($tpi,$tpj,$ln1,$ln2,@seq1,@seq2);
sub single {
        my ($seqf,$seqs)=@_;
        &cal($seqf,$seqs);
        my ($i,$j,$ii,$jj,$tt_max)=&trace_num;
        return ($i,$j,$ii,$jj,$tt_max);
}

sub twice {
        my ($seqf,$seqs)=@_;
        &cal($seqf,$seqs);
        my ($align1,$align2)=&trace_seq;
        return ($align1,$align2);
}

sub poscore{
        my ($seqf,$seqs)=@_;
        &cal($seqf,$seqs);
        my ($i,$j,$ii,$jj,$tt_max)=&trace_num;
        my $ln = 2*$ln2;
        return ($i, $tt_max/$ln);
}

sub score{
    my ($seqf,$seqs)=@_;
    &cal($seqf,$seqs);
    my $ln = $ln1 > $ln2 ? $ln2 : $ln1;
    $ln*=2;
    return $t_max/$ln;
}


sub cal{
        $t_max =0;
        my @tmp_arr = qw();
        @seq1=@tmp_arr; @seq2 =@tmp_arr; @Sco=@tmp_arr; @Poi=@tmp_arr; @Poj=@tmp_arr;
        my ($sq1,$sq2) = @_;
        $ln1=length $sq1;
        $ln2=length $sq2;
        @seq1 = split /\s*/,$sq1;
        @seq2 = split /\s*/,$sq2;
### intial the score and derived position matrix ###
        for my $i (0..$ln1){
                for my $j (0..$ln2){
                        $Sco[$i][$j]=0;
                        $Poi[$i][$j]=0;
                        $Poj[$i][$j]=0;
                }
        }

### initial the first row and first col ###
        for (my $i=1;$i<=$ln1;$i++){
                $Sco[$i][0]+=$Gap;
        }
        for my $j (1..$ln2){
                $Sco[0][$j]+=$Gap;
        }

### find the best score and record their path ###

        for my $i (1..$ln1){
                for my $j (1..$ln2){
                        my $m = ($seq1[$i-1] eq $seq2[$j-1])?1:0;
                        my ($ds,$di,$dj)= &mydia ($i,$j,$m);
                        my ($gs,$gi,$gj)= &mygap ($i,$j);
                        my ($is,$ii,$ij)= &myins ($i,$j);
                        if ($ds >= $gs){
                                if ($is > $ds){
                                        ($Sco[$i][$j],$Poi[$i][$j],$Poj[$i][$j]) = ($is,$ii,$ij);
                                }else{
                                        ($Sco[$i][$j],$Poi[$i][$j],$Poj[$i][$j]) = ($ds,$di,$dj);
                                }
                        }else{
                                if ($gs > $ii){
                                        ($Sco[$i][$j],$Poi[$i][$j],$Poj[$i][$j]) = ($gs,$gi,$gj);
                                }else{
                                        ($Sco[$i][$j],$Poi[$i][$j],$Poj[$i][$j]) = ($is,$ii,$ij);
                                }
                        } 
                        ($t_max,$tpi,$tpj) = ($Sco[$i][$j],$i,$j) if ($Sco[$i][$j] >= $t_max);
                }
        }
}
sub trace_num{
        my ($mi,$mj)=($tpi,$tpj);
        my ($i,$j);
        TNN:while (1){
                $i=$Poi[$tpi][$tpj];
                $j=$Poj[$tpi][$tpj];
                $tpi=$i;
                $tpj=$j;
                last TNN if ($i==0 or $j==0);
        }
        return ($i,$mi,$j,$mj,$t_max);
}

sub trace_seq{
#### trace back the aligned sequences ####
        my $aln1;
        my $aln2;
        my $ie=$ln1-$tpi;
        my $je=$ln2-$tpj;
        if ($ie < $je){
                for (my $j=$ln2; $j > $tpj; $j--){
                        $aln2.=$seq2[$j-1];
                        my $j_i=$j-$tpj-$ie;
                        if ($j_i > 0){$aln1.="-"}else{$aln1.=$seq1[$j_i-1]}
                }
        }else{
                for (my $i=$ln1; $i > $tpi; $i--){
                        $aln1.=$seq1[$i-1];
                        my $j_j=$i-$tpi-$je;
                        if ($j_j > 0){$aln2.="-"}else{$aln2.=$seq2[$j_j-1]}
                }
        }
        TNN:while (1){
                my $i=$Poi[$tpi][$tpj];
                my $j=$Poj[$tpi][$tpj];
                my $b1= $seq1[$tpi-1];
                my $b2= $seq2[$tpj-1];
                if ($i < $tpi){
                        $aln1.=$b1;
                }else{
                        $aln1.="-";
                }
                if ($j < $tpj){
                        $aln2.=$b2;
                }else{
                        $aln2.="-";
                }
                $tpi=$i;
                $tpj=$j;
                last TNN if ($i==0 or $j==0);
        }
        if ($tpi==0){
                for (my $jj=$tpj; $jj>0; $jj--){
                        $aln2.=$seq2[$jj-1];
                         $aln1.="-";
                }
        }elsif($tpj==0){
                for (my $ii=$tpi; $ii>0;$ii--){
                        $aln1.=$seq1[$ii-1]; 
                        $aln2.="-";
                }
        }

        $aln1=reverse $aln1;
        $aln2=reverse $aln2;
        return ($aln1,$aln2);
}

sub mydia {
        my ($si,$sj,$match)=@_;
        my $new = ($match==1 ? ($Sco[$si-1][$sj-1] + $Mat) : ($Sco[$si-1][$sj-1] - $Mis));
        return ($new, $si-1, $sj-1);
}

sub mygap {
        my ($si,$sj)=@_;
        my $max=$Sco[$si][$sj-1]-$Gap;
        my $spj=0;
        for my $k (0..($sj-1)){
                my $n=$Sco[$si][$k]-$Gap-(($sj-$k-1)*$Gex);
                ($max,$spj)= $max>$n?($max,$k):($n,$k);
        }
        return ($max,$si,$spj);
}

sub myins {
        my ($si,$sj)=@_;
        my $max=$Sco[$si-1][$sj]-$Gap;
        my $spi=0;
        for my $k (0..($si-1)){
                my $n=$Sco[$k][$sj]-$Gap-(($si-$k-1)*$Gex);
                ($max,$spi)= $max>$n?($max,$k):($n,$k);
        }
        return ($max,$spi,$sj);  
}

1;

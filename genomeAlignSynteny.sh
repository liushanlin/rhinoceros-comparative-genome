#! /bin/bash
if [[ $# -ne 7 ]]; then
	echo "usage bash genomeAlign.sh mafFile targeted.size query.size targeted.2bit query.2bit targetedPrefix queryPrefix"
	exit 0
fi
in="$1"
pre=$(echo $in | perl -pe 'chomp; s/\.maf//')
psl=${pre}.psl
ts="$2"
qs="$3"
tb="$4"
qb="$5"
tp="$6"
qp="$7"
maf-convert psl $in >$psl
chain=${pre}.chain.gz
axtChain -psl -linearGap=loose $psl $tb $qb stdout |gzip -c > $chain
net=${pre}.net.gz
chainPreNet $chain $ts $qs stdout | chainNet -minSpace=1 stdin $ts $qs stdout /dev/null | netSyntenic stdin stdout |gzip -c > $net
synNet=${pre}_syn.net.gz
netFilter -syn $net |gzip -c > $synNet
out=${pre}_syn.maf.gz
netToAxt $synNet $chain $tb $qb stdout |axtSort stdin stdout |axtToMaf -tPrefix=${tp}. -qPrefix=${qp}. stdin $ts $qs stdout |gzip -c > $out




#!/bin/zsh

# Handy shortcuts for working with antibodies, antibody repertoires, and cluster job evaluation

#alias calcdecoytime='grep attempted out | awk "{j+=\$10;n+=\$6;print n,\$10, j, j/n;}"'
function calcdecoytime(){
	grep attempted $@ | awk "{j+=\$10;n+=\$6;print n,\$10, j, j/n;}"
}

# note: ${1+$1/} expands to $1/ if $1 exists, otherwise empty.  allows a directory prefix or not
function lastpdb() {
	for d in ${1+$1/}*/pdbs; do echo -n $d/; ls $d | tail -1; done	
}

function H3lengths() {
	for f in ${1+$1/}*/grafting/details/H3.fasta; do echo -n $f\ ;tail -1 $f | wc | awk {'print $3-1'} ; done
}

function H3list() {
	for f in ${1+$1/}*/grafting/details/H3.fasta; do echo -n $f\ ;tail -1 $f; done 
}

function H3lengthshist() {
	for f in ${1+$1/}*/grafting/details/H3.fasta; do tail -1 $f | wc | awk {'print $3-1'} ; done | sort -n | uniq -c
}

function L1lengths() {
	for f in ${1+$1/}*/grafting/details/L1.fasta; do echo -n $f\ ;tail -1 $f | wc | awk {'print $3-1'} ; done
}

function L1list() {
	for f in ${1+$1/}*/grafting/details/L1.fasta; do echo -n $f\ ;tail -1 $f; done 
}

function L1lengthshist() {
	for f in ${1+$1/}*/grafting/details/L1.fasta; do tail -1 $f | wc | awk {'print $3-1'} ; done | sort -n | uniq -c
}

function Abjobtime() {
	for d in ${1+$1/}*/outerr
	do
		dirn=`dirname $d`
		mins=$(grep real $d/*err | cut -d'm' -f1 | cut -f2 | awk '{sum+=$1} END {print sum}')
		timeout=$(cat $d/*err | grep -c -i 'time limit')
		echo $dirn: $mins minutes + $timeout timeouts
	done
	echo
	grep 'SBATCH -t' ${1+$1/}abH3.sbatch
}

#!/bin/bash
if [ $# -ne 2 ] && [ $# -ne 3 ];then
	echo "Argument: breakdancer.SV.output size_low(bp) [size_up(<=3e9)]" && exit 1
fi 

bd=$1
size_low=$2
size_up=3000000000
size_suff="size_${size_low}"
if [ $# -eq 3 ];then
	size_up=$3
	size_suff="size_${size_low}_${size_up}"
fi
Dir=$(dirname $1)
BN=$(basename $bd .SV.output)

awk '{
	if($0 ~ /^#/){
	print
	}else{
		if(($7 == "ITX") || ($7 == "CTX")){
			if(($5-$2 >= LOW) && ($5-$2 <= UP)){
				print
			}
		}else if(($7 == "INS") || ($7 == "INV") || ($7 == "DEL")){
			SIZE=$8
			if(SIZE < 0){SIZE = (-1)*SIZE}
			if((SIZE >= LOW) && (SIZE <= UP)){
				print
			}
		}else{
			print "UNKNOWN TYPE"
		}
	}
}' LOW=$size_low UP=$size_up $bd > "${Dir}/${BN}.${size_suff}.SV.output"

if grep "UNKNOWN TYPE" "${Dir}/${BN}.${size_suff}.SV.output" > /dev/null;then 
	echo "Check ${Dir}/${BN}.${size_suff}.SV.output by 'UNKNOWN TYPE'"
fi


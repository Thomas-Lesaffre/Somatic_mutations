#!/bin/sh

if ls SAVE_* 1>/dev/null 2>&1; then
	ls SAVE_* > tmp
	s=$(cat tmp | sed 's|[^0-9]||g' | sort -n | uniq | tail -n 1)
else
	s=0
fi

cp s_parameters.txt parameters.txt
sed -i 's/SAUVEGARDE/'${s}'/g' parameters.txt
rm tmp	

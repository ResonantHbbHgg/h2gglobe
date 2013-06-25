#!/bin/bash
# to produce file list as input to lumi calculation
# O. Bondu, April 2013

dirFile="dirList.txt"
fileList="fileList.txt"

if [[ -e ${fileList} ]]
then
	echo "removing ${fileList}"
	rm ${fileList}
fi

for dir in `cat ${dirFile}`
do
	echo "listing files in ${dir}"
	for file in `cmsLs ${dir} | awk '{print $5}'`
	do
		cmsPfn ${file} >> ${fileList}
	done
done

exit 0

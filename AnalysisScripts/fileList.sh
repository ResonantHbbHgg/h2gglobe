#!/bin/bash
# to produce file list as input to lumi calculation
# O. Bondu, April 2013

eoscommand="/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select"

dirFile=""
if [[ -z ${1} ]]
then
	dirFile="dirList.txt"
fi

dirFile="${1}"


fileList="fileList.txt"

if [[ -e ${fileList} ]]
then
	echo "removing ${fileList}"
	rm ${fileList}
fi

for dir in `cat ${dirFile}`
do
	echo "listing files in ${dir}"
#	for file in `cmsLs ${dir} | awk '{print $5}'`
	for file in `${eoscommand} ls ${dir}`
	do
#		cmsPfn ${file} >> ${fileList}
		echo "root://eoscms//eos/cms${dir}/${file}" >> ${fileList}
	done
done

#exit 0

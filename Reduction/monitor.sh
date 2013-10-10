#!/bin/bash

#./submit_all_radion.sh

for iter in `seq 1 120`
do
	echo "iteration ${iter}"
	date
	totjobs=`bjobs | wc -l`
	runjobs=`bjobs -r | wc -l`
	run1nh=`bjobs -r | grep 1nh | wc -l`
	run8nh=`bjobs -r | grep 8nh | wc -l`
	run1nd=`bjobs -r | grep 1nd | wc -l`
	penjobs=`bjobs -p | wc -l`
	echo "Currently ${totjobs} jobs on batch (${runjobs} running [${run1nh} on 1nh, ${run8nh} on 8nh, ${run1nd} on 1nd] ; ${penjobs} pending)"
	echo "Checking status + resubmission (if any)"
	for batch in `echo "mc_radion_signal data2012_RERECO mc_Summer12_RD1 mc_radion_SM-HH"`
	do
		./check_status.py ${batch} > ${batch}_status.txt
		cat ${batch}_status.txt | grep "submit_reduction.sh" > to_resubmit.sh
		cat to_resubmit.sh
		bash to_resubmit.sh
	done
#	grep "submit_reduction.sh" *v2_status.txt
	echo "Now sleep for 5 minutes"
	sleep 300
done


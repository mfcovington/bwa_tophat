#!/usr/bin/env sh
# bwa_tophat_cleanup.sh
# Mike Covington
# created: 2012-02-12
#
# Description: Removes large, unnecessary and easy-to-regenerate files created by bwa_tophat.pl
# USAGE: bwa_tophat_cleanup.sh /path/to/ouput/directory/
#
if [ $1 ] && [ -d $1 ]; then
	out_dir=$1
	for IL in `ls $out_dir`
	do
		if [ -f $out_dir/$IL/bwa/bwa.*.sai ]; then
			[ -f $out_dir/$IL/tophat/accepted_hits.bam ] && [ -f $out_dir/$IL/merged/bwa_tophat_*.sorted.dupl_rm.bam.bai ] &&
			rm $out_dir/$IL/bwa/bwa.*.sai &&
			rm $out_dir/$IL/bwa/bwa.*.sam &&
			rm $out_dir/$IL/bwa/bwa.*.sorted.bam &&
			rm $out_dir/$IL/bwa/bwa.*.sorted.dupl_rm.bam &&
			rm $out_dir/$IL/bwa/unmapped.bwa.*.sorted.dupl_rm.fq &&
			rm $out_dir/$IL/merged/bwa_tophat_*.sorted.bam &&
			echo -e "$IL:\tUnnecessary files removed" ||
			echo -e "$IL:\tSomething is wrong"
		else
			echo -e "$IL:\tUnnecessary files have already been removed"
		fi
	done
else
	echo ""
	echo "You did not specify a directory that exists."
	echo "USAGE: bwa_tophat_cleanup.sh /path/to/ouput/directory/"
	echo ""
fi
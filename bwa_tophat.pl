#!/usr/bin/perl
# bwa_tophat.pl
# Mike Covington
# created: 2012-02-04
#
# Description: 
#
use strict; use warnings;
use Getopt::Long;

my ($fq_id, @fq_in, $fq_cat, $multi_fq, $out_dir);

### GETOPT DEFAULTS
my $threads = 1;
my $bwa_e = 15;
my $bwa_i = 10;
my $bwa_k = 1;
my $bwa_l = 25;
my $bwa_n = 0.02;
my $bwa_db = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/BWA_S_lycopersicum_chromosomes.2.40";
my $ref_id = "Slyc";
my $samse_n = 1;
my $bowtie_db = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/Bowtie_S_lycopersicum_chromosomes.2.40";
my $ref_fa = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa";



my $options = GetOptions (
	"out=s"		=>	\$out_dir,
	"threads=i"	=>	\$threads,
	"bwa_e=i"	=>	\$bwa_e,
	"bwa_i=i"	=>	\$bwa_i,
	"bwa_k=i"	=>	\$bwa_k,
	"bwa_l=i"	=>	\$bwa_l,
	"bwa_n=f"	=>	\$bwa_n,
	"bwa_db=s"	=>	\$bwa_db,
	"fq=s{1,}"	=>	\@fq_in,
	"fq_id=s"	=>	\$fq_id,
	"ref_id=s"	=>	\$ref_id,
	"samse_n=i"	=>	\$samse_n,
	"bowtie_db=s"	=>	\$bowtie_db,
	"ref_fa=s"	=>	\$ref_fa
);


my $log_dir = "$out_dir/logs/";
my $bwa_dir = "$out_dir/bwa/";
my $tophat_dir = "$out_dir/tophat/";
my $merged_dir = "$out_dir/merged/";

`mkdir -p $out_dir` unless -e $out_dir;
`mkdir $log_dir` unless -e $log_dir;
`mkdir $bwa_dir` unless -e $bwa_dir;
`mkdir $tophat_dir` unless -e $tophat_dir;
`mkdir $merged_dir` unless -e $merged_dir;


print "BEGIN!\n\n";

if (scalar @fq_in > 1) {
	$multi_fq = 1;
	print `date`;
	print"CMD: cat @fq_in > $out_dir/temp.fq\n\n";
	`cat @fq_in > $out_dir/temp.fq`;
	$fq_cat = "$out_dir/temp.fq";
}else{
	$fq_cat = $fq_in[0];
}

# align with bwa
print `date`;
print "CMD: bwa aln -t $threads -e $bwa_e -i $bwa_i -k $bwa_k -l $bwa_l -n $bwa_n $bwa_db $fq_cat > $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sai 2> $log_dir/$fq_id-$ref_id.bwa.aln.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.log\n\n";
`bwa aln -t $threads -e $bwa_e -i $bwa_i -k $bwa_k -l $bwa_l -n $bwa_n $bwa_db $fq_cat > $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sai 2> $log_dir/$fq_id-$ref_id.bwa.aln.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.log`;

#convert sai to sam
print `date`;
print "CMD: bwa samse -n $samse_n $bwa_db $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sai $fq_cat > $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sam 2> $log_dir/$fq_id-$ref_id.bwa.samse.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.log\n\n";
`bwa samse -n $samse_n $bwa_db $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sai $fq_cat > $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sam 2> $log_dir/$fq_id-$ref_id.bwa.samse.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.log`;

#convert sam to bam
print `date`;
print "CMD: samtools view -Shb -o $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.bam $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sam 2> $log_dir/bwa_sam2bam.log\n\n";
`samtools view -Shb -o $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.bam $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sam 2> $log_dir/$fq_id-$ref_id.bwa_sam2bam.log`;

####sort
print `date`;
print "CMD: java -jar /usr/local/bin/SortSam.jar INPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.bam OUTPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> $log_dir/$fq_id-$ref_id.sort.1.log\n\n";
`java -jar /usr/local/bin/SortSam.jar INPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.bam OUTPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> $log_dir/$fq_id-$ref_id.sort.1.log`;

####remove duplicates
print `date`;
print "CMD: java -jar /usr/local/bin/MarkDuplicates.jar INPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.bam OUTPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam METRICS_FILE=$bwa_dir/$fq_id-$ref_id.picard.dupl_rm.metrics_file.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $log_dir/$fq_id-$ref_id.rm_dup.1.log\n\n";
`java -jar /usr/local/bin/MarkDuplicates.jar INPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.bam OUTPUT=$bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam METRICS_FILE=$bwa_dir/$fq_id-$ref_id.picard.dupl_rm.metrics_file.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $log_dir/$fq_id-$ref_id.rm_dup.1.log`;

####extract mapped reads to bam file
print `date`;
print "CMD: samtools view -hb -F 4 $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam > $bwa_dir/mapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam 2> $log_dir/$fq_id-$ref_id.extract_bwa_mapped.log\n\n";
`samtools view -hb -F 4 $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam > $bwa_dir/mapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam 2> $log_dir/$fq_id-$ref_id.extract_bwa_mapped.log`;

#extract unmapped reads to bam file
print `date`;
print "CMD: samtools view -hb -f 4 $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam > $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam 2> $log_dir/$fq_id-$ref_id.extract_bwa_unmapped.log\n\n";
`samtools view -hb -f 4 $bwa_dir/bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam > $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam 2> $log_dir/$fq_id-$ref_id.extract_bwa_unmapped.log`;

#convert unmapped.bam to unmapped.fq
print `date`;
print "CMD: samtools bam2fq $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam > $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.fq 2> $log_dir/$fq_id-$ref_id.bam2fq.log\n\n";
`samtools bam2fq $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam > $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.fq 2> $log_dir/$fq_id-$ref_id.bam2fq.log`;

#align unmapped with tophat
print `date`;
print "CMD: tophat --splice-mismatches 1 --max-multihits 1 --segment-length 22 --butterfly-search --max-intron-length 5000 --solexa1.3-quals --num-threads $threads --library-type fr-unstranded --output-dir $tophat_dir $bowtie_db $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.fq 2> $log_dir/tophat.$fq_id-$ref_id.tophat_unmapped.sm_1.mm_1.sl_22.mil_5000.lt_fu.bs.s13q.log\n\n";
`tophat --splice-mismatches 1 --max-multihits 1 --segment-length 22 --butterfly-search --max-intron-length 5000 --solexa1.3-quals --num-threads $threads --library-type fr-unstranded --output-dir $tophat_dir $bowtie_db $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.fq 2> $log_dir/$fq_id-$ref_id.tophat_unmapped.sm_1.mm_1.sl_22.mil_5000.lt_fu.bs.s13q.log`;

#extract reads still unmapped after tophat
#convert to sam
print `date`;
print "CMD: Extracting unmapped tophat reads, etc.\n\n";
`samtools view -h $tophat_dir/accepted_hits.bam > $tophat_dir/accepted_hits.sam;   cp -p $tophat_dir/accepted_hits.sam $tophat_dir/accepted_hits_and_unmapped.sam;   /Volumes/OuterColomaData/Mike/extract_unmapped_tophat_to_sam.pl --fq $bwa_dir/unmapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.fq --sam $tophat_dir/accepted_hits.sam --out ">$tophat_dir/accepted_hits_and_unmapped.sam";   samtools view -Shb $tophat_dir/accepted_hits_and_unmapped.sam > $tophat_dir/accepted_hits_and_unmapped.bam;   rm $tophat_dir/accepted_hits.sam $tophat_dir/accepted_hits_and_unmapped.sam`;

#####changed this merge to merge bwa_mapped + tophat_mapped + tophat_unmapped
#merge bwa and tophat bam files
print `date`;
print "CMD: samtools merge $merged_dir/bwa_tophat_$fq_id-$ref_id.bam $bwa_dir/mapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam $tophat_dir/accepted_hits_and_unmapped.bam > $log_dir/$fq_id-$ref_id.merge.log\n\n";
`samtools merge $merged_dir/bwa_tophat_$fq_id-$ref_id.bam $bwa_dir/mapped.bwa.$fq_id-$ref_id.e_$bwa_e.i_$bwa_i.k_$bwa_k.l_$bwa_l.n_$bwa_n.sorted.dupl_rm.bam $tophat_dir/accepted_hits_and_unmapped.bam > $log_dir/$fq_id-$ref_id.merge.log`;

#sort
print `date`;
print "CMD: java -jar /usr/local/bin/SortSam.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> $log_dir/$fq_id-$ref_id.sort.2.log\n\n";
`java -jar /usr/local/bin/SortSam.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> $log_dir/$fq_id-$ref_id.sort.2.log`;

#remove duplicates  (can I remove this since I rm_dup'd earlier?)
print `date`;
print "CMD: java -jar /usr/local/bin/MarkDuplicates.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam METRICS_FILE=$merged_dir/$fq_id-$ref_id.picard.dupl_rm.metrics_file.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $log_dir/$fq_id-$ref_id.rm_dup.2.log\n\n";
`java -jar /usr/local/bin/MarkDuplicates.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam METRICS_FILE=$merged_dir/$fq_id-$ref_id.picard.dupl_rm.metrics_file.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $log_dir/$fq_id-$ref_id.rm_dup.2.log`;

#index
print `date`;
print "CMD: java -jar /usr/local/bin/BuildBamIndex.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam.bai 2> $log_dir/$fq_id-$ref_id.index.log\n\n";
`java -jar /usr/local/bin/BuildBamIndex.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam.bai 2> $log_dir/$fq_id-$ref_id.index.log`;

#index
# print "\nCMD: java -jar /usr/local/bin/BuildBamIndex.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam.bai 2> $log_dir/$fq_id-$ref_id.index.2.log\n\n";
# `java -jar /usr/local/bin/BuildBamIndex.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam.bai 2> $log_dir/$fq_id-$ref_id.index.2.log`;

# print "\nCMD: java -jar /usr/local/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref_fa -o $merged_dir/$fq_id-$ref_id.target_positions.list -I $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam  > $log_dir/$fq_id-$ref_id.create_target.log\n\n";
# `java -jar /usr/local/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref_fa -o $merged_dir/$fq_id-$ref_id.target_positions.list -I $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam  > $log_dir/$fq_id-$ref_id.create_target.log`;

# print "\nCMD: java -Djava.io.tmpdir=$out_dir -jar /usr/local/bin/GenomeAnalysisTK.jar -I $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam -R $ref_fa -T IndelRealigner -targetIntervals $merged_dir/$fq_id-$ref_id.target_positions.list -o $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.bam > $log_dir/$fq_id-$ref_id.realign.log\n\n";
# `java -Djava.io.tmpdir=$out_dir -jar /usr/local/bin/GenomeAnalysisTK.jar -I $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam -R $ref_fa -T IndelRealigner -targetIntervals $merged_dir/$fq_id-$ref_id.target_positions.list -o $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.bam > $log_dir/$fq_id-$ref_id.realign.log`;

# print "\nCMD: java -jar /usr/local/bin/SortSam.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.bam OUTPUT= $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> $log_dir/$fq_id-$ref_id.sort.2.log\n\n";
# `java -jar /usr/local/bin/SortSam.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.bam OUTPUT= $merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> $log_dir/$fq_id-$ref_id.sort.2.log`;

# print "\nCMD: java -jar /usr/local/bin/BuildBamIndex.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.sorted.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.sorted.bam.bai 2> $log_dir/$fq_id-$ref_id.index.2.log\n\n";
# `java -jar /usr/local/bin/BuildBamIndex.jar INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.sorted.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.realigned.sorted.bam.bai 2> $log_dir/$fq_id-$ref_id.index.2.log`;

`rm $out_dir/temp.fq` if $multi_fq ==1;

print "END!\n";
`date`;

exit;













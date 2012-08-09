#!/usr/bin/perl
# bwa_tophat.pl
# Mike Covington
# created: 2012-02-04
#
# Description: 
#
use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use File::Path 'make_path';
use File::Copy 'cp';
use IPC::Cmd qw[can_run run run_forked];
use DateTime;

use bwa_aln_commander

#TODO: get rid of concatenating fq files!
#TODO: give fa instead of prefix for bowtie/bwa index; check for existence of index files and create, if necessary
#TODO: for bowtie (and others?) have an extra argument for parameters not supported by my module

my ($fq_id, @fq_in, $multi_fq, $out_dir, $verbose);

### GETOPT DEFAULTS
my $threads = 1;
my $bwa_e = 15;
my $bwa_i = 10;
my $bwa_k = 1;
my $bwa_l = 25;
my $bwa_n = 0.05;
my $bwa_db = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/BWA_S_lycopersicum_chromosomes.2.40";
my $ref_id = "Slyc";
my $samse_n = 1;
my $bowtie_db = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/Bowtie_S_lycopersicum_chromosomes.2.40";
my $ref_fa = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa";



my $options = GetOptions(
    "out=s"       => \$out_dir,
    "threads=i"   => \$threads,
    "bwa_e=i"     => \$bwa_e,
    "bwa_i=i"     => \$bwa_i,
    "bwa_k=i"     => \$bwa_k,
    "bwa_l=i"     => \$bwa_l,
    "bwa_n=f"     => \$bwa_n,
    "bwa_db=s"    => \$bwa_db,
    "fq=s{1,}"    => \@fq_in,
    "fq_id=s"     => \$fq_id,
    "ref_id=s"    => \$ref_id,
    "samse_n=i"   => \$samse_n,
    "bowtie_db=s" => \$bowtie_db,
    "ref_fa=s"    => \$ref_fa,
    "verbose"     => \$verbose,
);

#check for dependencies
my @bins = qw( bwa samtools tophat java );
my @jars  = qw( /usr/local/bin/SortSam.jar /usr/local/bin/MarkDuplicates.jar /usr/local/bin/BuildBamIndex.jar );
can_run($_) or die "$_ is not installed." for @bins;
-e $_       or die "$_ is not installed." for @jars;
my ( $cmd, $tool );

#set parameters
my $bwa_aln_param = "-t $threads -e $bwa_e -i $bwa_i -k $bwa_k -l $bwa_l -n $bwa_n $bwa_db $fq_cat > $bwa_dir/bwa.$fq_id-$ref_id.sai";
my $bwa_samse_param = "-n $samse_n $bwa_db $bwa_dir/bwa.$fq_id-$ref_id.sai $fq_cat > $bwa_dir/bwa.$fq_id-$ref_id.sam";
my $sam2bam_param_1 = "-Shb -o $bwa_dir/bwa.$fq_id-$ref_id.bam $bwa_dir/bwa.$fq_id-$ref_id.sam";
my $sortsam_param_1 = "INPUT=$bwa_dir/bwa.$fq_id-$ref_id.bam OUTPUT=$bwa_dir/bwa.$fq_id-$ref_id.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT";
my $markdups_param_1 = "INPUT=$bwa_dir/bwa.$fq_id-$ref_id.sorted.bam OUTPUT=$bwa_dir/bwa.$fq_id-$ref_id.sorted.dupl_rm.bam METRICS_FILE=$bwa_dir/$fq_id-$ref_id.picard.dupl_rm.metrics_file.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true";
my $extract_mapped_param = "-hb -F 4 $bwa_dir/bwa.$fq_id-$ref_id.sorted.dupl_rm.bam > $bwa_dir/mapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.bam";
my $extract_unmapped_param = "-hb -f 4 $bwa_dir/bwa.$fq_id-$ref_id.sorted.dupl_rm.bam > $bwa_dir/unmapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.bam";
my $bam2fq_param = "$bwa_dir/unmapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.bam > $bwa_dir/unmapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.fq";
my $tophat_param = "--splice-mismatches 1 --max-multihits 1 --segment-length 22 --butterfly-search --max-intron-length 5000 --solexa1.3-quals --num-threads $threads --library-type fr-unstranded --output-dir $tophat_dir $bowtie_db $bwa_dir/unmapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.fq";

my $sam2bam_param_2 = "-Shb $tophat_dir/accepted_hits_and_unmapped.sam > $tophat_dir/accepted_hits_and_unmapped.bam";


my $merge_param = "$merged_dir/bwa_tophat_$fq_id-$ref_id.bam $bwa_dir/mapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.bam $tophat_dir/accepted_hits_and_unmapped.bam";
my $sortsam_param_2 = "INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT";
my $markdups_param_2 = "INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam METRICS_FILE=$merged_dir/$fq_id-$ref_id.picard.dupl_rm.metrics_file.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true";
my $bam_index_param = "INPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam OUTPUT=$merged_dir/bwa_tophat_$fq_id-$ref_id.sorted.dupl_rm.bam.bai";


my $bwa_aln = bwa_aln_commander->new(
    in_fq     => $fq_cat,
    ref_fasta => $ref_fasta,
    threads   => $threads,
    e         => $bwa_e,
    i         => $bwa_i,
    k         => $bwa_k,
    l         => $bwa_l,
    n         => $bwa_n,
    file_out  => "$bwa_dir/bwa.$fq_id-$ref_id.sai",
);



#make output directories
my $log_dir = "$out_dir/logs/";
my $bwa_dir = "$out_dir/bwa/";
my $tophat_dir = "$out_dir/tophat/";
my $merged_dir = "$out_dir/merged/";
make_path( $out_dir, $log_dir, $bwa_dir, $tophat_dir, $merged_dir );

open my $stdout_log, ">", "$log_dir/out.log";
open my $stderr_log, ">", "$log_dir/err.log";
open my $cmd_log,    ">", "$log_dir/cmd.log";

my $fq_cat;
if (scalar @fq_in > 1) {
    $multi_fq = 1;
    $fq_cat = "$out_dir/temp.fq";
    run_cmd( [ "cat", @fq_in, ">", $fq_cat ] );
}else{
    $fq_cat = $fq_in[0];
}

run_cmd( [ "bwa aln", $bwa_aln_param ] );
run_cmd( [ "bwa samse", $bwa_samse_param ] );

# filter sam for -q30 and XT:A:U


# how to keep track of reads that get filtered out so they can be sent to tophat??

# grep -ve "XT:A:U"

run_cmd( [ "samtools view",                               $sam2bam_param_1          ] );
run_cmd( [ "java -jar /usr/local/bin/SortSam.jar",        $sortsam_param_1        ] );
run_cmd( [ "java -jar /usr/local/bin/MarkDuplicates.jar", $markdups_param_1       ] );
run_cmd( [ "samtools view",                               $extract_mapped_param   ] );
run_cmd( [ "samtools view",                               $extract_unmapped_param ] );
run_cmd( [ "samtools bam2fq",                             $bam2fq_param           ] );
run_cmd( [ "tophat",                                      $tophat_param           ] );

#extract reads still unmapped after tophat
#convert to sam
print `date`;
print "CMD: Extracting unmapped tophat reads, etc.\n\n";
`samtools view -h $tophat_dir/accepted_hits.bam > $tophat_dir/accepted_hits.sam`;

cp "$tophat_dir/accepted_hits.sam", "$tophat_dir/accepted_hits_and_unmapped.sam" or die "Can't copy";
`/Volumes/OuterColomaData/Mike/extract_unmapped_tophat_to_sam.pl --fq $bwa_dir/unmapped.bwa.$fq_id-$ref_id.sorted.dupl_rm.fq --sam $tophat_dir/accepted_hits.sam --out ">$tophat_dir/accepted_hits_and_unmapped.sam"`;
run_cmd( [ "samtools view",                               $sam2bam_param_2          ] );
unlink "$tophat_dir/accepted_hits.sam", "$tophat_dir/accepted_hits_and_unmapped.sam";

#####changed this merge to merge bwa_mapped + tophat_mapped + tophat_unmapped
#merge bwa and tophat bam files
run_cmd( [ "samtools merge",                               $merge_param ] );
run_cmd( [ "java -jar /usr/local/bin/SortSam.jar",        $sortsam_param_2        ] );
run_cmd( [ "java -jar /usr/local/bin/MarkDuplicates.jar", $markdups_param_2       ] );
run_cmd( [ "java -jar /usr/local/bin/BuildBamIndex.jar", $bam_index_param       ] );
unlink "rm $out_dir/temp.fq" if $multi_fq ==1;

sub run_cmd {
    my $cmd = shift;
    my $datestamp = DateTime->now . " ---- " . join " ", @$cmd;
    say $stdout_log $datestamp;
    say $stderr_log $datestamp;
    say $cmd_log $datestamp;
    my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $cmd, verbose => $verbose );
    die $error_message unless $success;
    print $stdout_log join "", @$stdout_buf;
    print $stderr_log join "", @$stderr_buf;
}

exit;













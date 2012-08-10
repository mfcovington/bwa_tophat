#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2012-08-09
#
# Description: 
#
use strict;
use warnings;
use autodie;
use feature 'say';
use bwa_aln_commander;
use Data::Printer;

my $bwa = bwa_aln_commander->new(
    fq_id => "WT",
    ref_id => "organism",
    fq_in => "sequences.fq",
    fasta_ref => "/refs/reference.fasta",
    n => 0.04,
    o => 55,
    se => 1,
    pe2 => 1,
    verbose => 1,
    # _log_dir => "kellogs2",
    # out_dir => "",
    # _stdout_log => "test_stdout_log",
);

say $bwa->_param;
say $bwa->_tool;
say $bwa->_cmd;
p $bwa->_cmd;

$bwa->_run_cmd;

# $bwa->_tool('echo');
# say $bwa->_tool;
# say 'hi';
# $bwa->_run_cmd;

# say "hi";
# $bwa->_open_fhs;
# say $bwa->_cmd;


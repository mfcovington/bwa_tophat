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
    fq_in => "sequences.fq",
    fasta_ref => "/refs/reference.fasta",
    n => 0.04,
    o => 55,
    se => 1,
    pe2 => 1,
    verbose => 1,
    out_dir => "~/sandbox",
    _stdout_log => "test_stdout_log",
);

say $bwa->_param;
say $bwa->_tool;
say $bwa->_cmd;
p $bwa->_cmd;

$bwa->_run_cmd;
say "hi";
# $bwa->_open_fhs;
# say $bwa->_cmd;


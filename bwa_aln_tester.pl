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

my $bwa = bwa_aln_commander->new(
    fq_in => "sequences.fq",
    n => 0.04,
    o => 55,
);

say $bwa->_bwa_aln_param;

package bwa_aln_commander;
use Moose;
use MooseX::Aliases;
use MooseX::Clone;
use MooseX::UndefTolerant;
use File::Basename;
# use File::Path 'make_path';
use Parallel::ForkManager;
use DateTime;
# use IPC::Cmd qw[can_run run run_forked];
use autodie;
use feature 'say';
with qw(MooseX::Clone);
extends 'commander';

use Data::Printer;

###TODO: make auto-generate f from fq_in? (only if not specified?)
###TOTO: at some point make a check for required params: prefix and fq_in
###TODO: prevent -0, -1, and -2 together? (AND make sure -b is used, too?)
###TODO: validity tests can_run bwa, test version, test flag compatibilities (like -012b mentionsabove)
###TODO: make sure actually using all modules loaded

sub _param {
    my $self = shift;

    my @args = ( 'n', 'o', 'e', 'i', 'd', 'l', 'k', 'm', 't', 'M', 'O', 'E', 'R', 'q', 'f', 'B', 'c', 'L', 'N', 'I', 'b' );
    my @params;
    for (@args) {
        push @params, join " ", "-$_", $self->$_ if defined $self->$_;
    }
    push @params, "-0" if $self->se;
    push @params, "-1" if $self->pe1;
    push @params, "-2" if $self->pe2;
    push @params, $self->_prefix;
    push @params, $self->fq;

    my $param_string = join " ", @params;
    return $param_string;
}

sub _prefix {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->fasta_ref, ".fa(sta)?" );
    return $dir_name . $filename;
}

# Don't clone these attributes
has 'fq_id' => (
    is    => 'rw',
    isa   => 'Str',
    traits => [qw(NoClone)],
);

has 'fq' => (
    is    => 'rw',
    isa   => 'Str',
    alias   => [qw(fq_in in_fq fastq)],
    traits => [qw(NoClone)],
    predicate => 'has_fq',
);

has 'f' => (
    is      => 'rw',
    isa     => 'Str',
    alias   => 'file_out',
    traits => [qw(NoClone)],
    lazy => 1,
    default => sub {
        my $self = shift;

        return $self->_bwa_dir . join ".", "bwa", $self->fq_id, $self->ref_id, "sai"
    },
    # -f FILE   file to write output to instead of stdout
);

has '_bwa_dir' => (
    is      => 'rw',
    isa     => 'Str',
    traits  => [qw(NoClone)],
    lazy => 1,
    default => sub {
        my $self = shift;

        return $self->out_dir . "/bwa/";
    },
);

# OK to clone these attributes
has '_tool' => (
    is    => 'ro',
    isa   => 'Str',
    default => 'bwa aln',
);

has 'fasta_ref' => (
    is    => 'rw',
    isa   => 'Str',
    alias   => 'ref_fasta',
);

has 'ref_id' => (
    is    => 'rw',
    isa   => 'Str',
);

has 'n' => (
    is    => 'rw',
    isa   => 'Num',
    alias => [qw(max_diff missing_prob)],
    # -n NUM    max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
);

has 'o' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'gap_open',
    # -o INT    maximum number or fraction of gap opens [1]
);

has 'e' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'gap_extend',
    # -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
);

has 'i' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'indel_distance',
    # -i INT    do not put an indel within INT bp towards the ends [5]
);

has 'd' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'max_deletion',
    # -d INT    maximum occurrences for extending a long deletion [10]
);

has 'l' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'seed_length',
    # -l INT    seed length [32]
);

has 'k' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'max_seed_diff',
    # -k INT    maximum differences in the seed [2]
);

has 'm' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'max_entries',
    # -m INT    maximum entries in the queue [2000000]
);

has 't' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'threads',
    # -t INT    number of threads [1]
);

has 'M' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'mismatch_penalty',
    # -M INT    mismatch penalty [3]
);

has 'O' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'gap_open_penalty',
    # -O INT    gap open penalty [11]
);

has 'E' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'gap_extend_penalty',
    # -E INT    gap extension penalty [4]
);

has 'R' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'max_hits',
    # -R INT    stop searching when there are >INT equally best hits [30]
);

has 'q' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'qual_cutoff',
    # -q INT    quality threshold for read trimming down to 35bp [0]
);

has 'B' => (
    is      => 'rw',
    isa     => 'Int',
    alias   => 'barcode_length',
    # -B INT    length of barcode
);

has 'c' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'color_space',
    # -c        input sequences are in the color space
);

has 'L' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'long_del_gap_penalty',
    # -L        log-scaled gap penalty for long deletions
);

has 'N' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'non_iterative',
    # -N        non-iterative mode: search for all n-difference hits (slooow)
);

has 'I' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'illumina_1.3_plus',
    # -I        the input is in the Illumina 1.3+ FASTQ-like format
);

has 'b' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'BAM_in',
    traits => [qw(NoClone)],
    # -b        the input read file is in the BAM format
);

has 'se' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'single',
    # -0        use single-end reads only (effective with -b)
);

has 'pe1' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'pair1',
    # -1        use the 1st read in a pair (effective with -b)
);

has 'pe2' => (
    is      => 'rw',
    isa     => 'Bool',
    alias   => 'pair2',
    # -2        use the 2nd read in a pair (effective with -b)
);

no Moose;
__PACKAGE__->meta->make_immutable;
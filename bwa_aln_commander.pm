package bwa_aln_commander;
use Moose;
use MooseX::Aliases;
use MooseX::Clone;
use MooseX::UndefTolerant;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use DateTime;
use IPC::Cmd qw[can_run run run_forked];
use autodie;
use feature 'say';
with qw(MooseX::Clone);

use Data::Printer;

###TODO: make auto-generate f from fq_in? (only if not specified?)
###TOTO: at some point make a check for required params: prefix and fq_in
###TODO: prevent -0, -1, and -2 together? (AND make sure -b is used, too?)
###TODO: specify out_dir
###TODO: validity tests can_run bwa, test version, test flag compatibilities (like -012b mentionsabove)
###TODO: make sure using all modules loaded

sub _tool {
    my $self = shift;

    return "bwa aln";
}

sub _open_fhs {
    my $self = shift;

    open my $fh, ">", $self->_stdout_log;
    # p $self->_stdout_log_fh; 
    # p $fh;
    $self->_stdout_log_fh($fh);
    # p $self->_stdout_log_fh; 
    say {$self->_stdout_log_fh} "opening!";
    # my $fh2 = $self->_stdout_log_fh;
    # say $fh2 "hi!";

}

sub _param {
    my $self = shift;

    my @args = ( 'n', 'o', 'e', 'i', 'd', 'l', 'k', 'm', 't', 'M', 'O', 'E', 'R', 'q', 'B', 'c', 'L', 'N', 'I', 'b' );
    my @params;
    for (@args) {
        push @params, join " ", "-$_", $self->$_ if defined $self->$_;
    }
    push @params, "-0" if $self->se;
    push @params, "-1" if $self->pe1;
    push @params, "-2" if $self->pe2;
    push @params, $self->_prefix;
    push @params, $self->fq;

    my $bwa_aln_param = join " ", @params;
    return $bwa_aln_param;
}

sub _cmd {
    my $self = shift;

    return [ $self->_tool, $self->_param ];
}

sub _run_cmd {
    my $self = shift;

    $self->_open_fhs;
    say {$self->_stdout_log_fh} "running!";
    p $self->_stdout_log_fh;
    my $cmd = $self->_cmd;
    my $datestamp = DateTime->now . " ---- " . join " ", @$cmd;
    # say $stdout_log $datestamp;
    # say $stderr_log $datestamp;
    # say $cmd_log $datestamp;
    my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $cmd, verbose => $self->verbose );
    die $error_message unless $success;
    # print $stdout_log join "", @$stdout_buf;
    # print $stderr_log join "", @$stderr_buf;
}

sub _prefix {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->fasta_ref, ".fa(sta)?" );
    return $dir_name . $filename;
}

# Don't clone these attributes
has 'fq' => (
    is    => 'rw',
    isa   => 'Str',
    alias => [qw(fastq fq)],
    alias   => [qw(fq_in in_fq)],
    traits => [qw(NoClone)],
    predicate => 'has_fq',
);

has 'f' => (
    is      => 'rw',
    isa     => 'Str',
    alias   => 'file_out',
    traits => [qw(NoClone)],
    # -f FILE   file to write output to instead of stdout
);

has 'fasta_ref' => (
    is    => 'rw',
    isa   => 'Str',
    alias   => 'ref_fasta',
);

has 'out_dir' => (
    is => 'rw',
    isa => 'Str',
    traits => [qw(NoClone)],
    default => "./",
);

has '_stdout_log' => (
    is => 'rw',
    isa => 'Str',
    traits => [qw(NoClone)],
);

has '_stdout_log_fh' => (
    is => 'rw',
    isa => 'FileHandle',
    traits => [qw(NoClone)],
);

has '_stderr_log' => (
    is => 'rw',
    isa => 'FileHandle',
    traits => [qw(NoClone)],
);

has '_cmd_log' => (
    is => 'rw',
    isa => 'FileHandle',
    traits => [qw(NoClone)],
);



# OK to clone these attributes
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

has 'verbose' => (
    is      => 'rw',
    isa     => 'Bool',
);

no Moose;
__PACKAGE__->meta->make_immutable;
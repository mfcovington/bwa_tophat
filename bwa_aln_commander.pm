package bwa_aln_commander;
use Moose;
use MooseX::Aliases;
use MooseX::Clone;
use MooseX::UndefTolerant;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
use feature 'say';
with qw(MooseX::Clone);

use Data::Printer;

###TODO: make auto-generate f from fq_in? (only if not specified?)
###TOTO: at some point make a check for required params: prefix and fq_in
###TODO: prevent -0, -1, and -2 together?
###TODO: specify out_dir

# sub _run_cmd {
#     my $cmd = shift;
#     my $datestamp = DateTime->now . " ---- " . join " ", @$cmd;
#     say $stdout_log $datestamp;
#     say $stderr_log $datestamp;
#     say $cmd_log $datestamp;
#     my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
#       run( command => $cmd, verbose => $verbose );
#     die $error_message unless $success;
#     print $stdout_log join "", @$stdout_buf;
#     print $stderr_log join "", @$stderr_buf;
# }

sub _bwa_aln_param {
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
    p @args;
    p @params;

    my $bwa_aln_param = join " ", @params;
    return $bwa_aln_param;
}


# sub _bwa_aln_param {
#     my $self = shift;

#     my %params;
#     $params{'n'} = $self->n   if $self->has_n;
#     $params{'o'} = $self->o   if $self->has_o;
    # $params{'e'} = $self->e   if $self->has_e;
    # $params{'i'} = $self->i   if $self->has_i;
    # $params{'d'} = $self->d   if $self->has_d;
    # $params{'l'} = $self->l   if $self->has_l;
    # $params{'k'} = $self->k   if $self->has_k;
    # $params{'m'} = $self->m   if $self->has_m;
    # $params{'t'} = $self->t   if $self->has_t;
    # $params{'M'} = $self->M   if $self->has_M;
    # $params{'O'} = $self->O   if $self->has_O;
    # $params{'E'} = $self->E   if $self->has_E;
    # $params{'R'} = $self->R   if $self->has_R;
    # $params{'q'} = $self->q   if $self->has_q;
    # $params{'B'} = $self->B   if $self->has_B;
    # $params{'c'} = $self->c   if $self->has_c;
    # $params{'L'} = $self->L   if $self->has_L;
    # $params{'N'} = $self->N   if $self->has_N;
    # $params{'I'} = $self->I   if $self->has_I;
    # $params{'b'} = $self->b   if $self->has_b;
    # $params{'0'} = $self->se  if $self->has_se;
    # $params{'1'} = $self->pe1 if $self->has_pe1;
    # $params{'2'} = $self->pe2 if $self->has_pe2;

#     my @param_array;
#     for (keys %params) {
#         push @param_array, join " ", "-$_", $params{$_};
#     }
#     push @param_array, "-0" if $self->has_se;
#     push @param_array, "-1" if $self->has_pe1;
#     push @param_array, "-2" if $self->has_pe2;
#     push @param_array, $self->prefix;
#     push @param_array, $self->fq;

#     my $bwa_aln_param = join " ", @param_array;
#     return $bwa_aln_param;
# }

# sub bwa_aln {
#     my $self = shift;

#     _run_cmd( [ "bwa aln", $bwa_aln_param ] );

#     my $bwa_aln_cmd = "bwa"
# }



#     extends 'Moose::Meta::Attribute';

# sub dump {
#     my $self = shift;

    # iterate over all the attributes in $self
    # my %attributes = %{ $self->meta->get_attribute_map };
    # while (my ($name, $attribute) = each %attributes) {

    #     # print the label if available
    #     if ($attribute->isa('MyApp::Meta::Attribute::Labeled')
    #         && $attribute->has_label) {
    #             print $attribute->label;
    #     }
    #     # otherwise print the name
    #     else {
    #         print $name;
    #     }

    #     # print the attribute's value
    #     my $reader = $attribute->get_read_method;
    #     print ": " . $self->$reader . "\n";
    # }
# }

sub samtools_cmd_gaps {
    my $self = shift;

    my $samtools_cmd = "samtools mpileup" . $self->_region . $self->bam . " | cut -f1-2,4 > " . $self->out_file . ".cov_gaps";
    return $samtools_cmd;
}

sub _prefix {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->fasta_ref, ".fa(sta)?" );
    return $dir_name . $filename;
}


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
# __PACKAGE__->meta->make_immutable;
__PACKAGE__->meta->make_immutable(inline_constructor => 0);


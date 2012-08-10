package commander;
use Moose;
use File::Path 'make_path';
use IPC::Cmd qw[can_run run run_forked];
use MooseX::Clone;
use MooseX::UndefTolerant;
use autodie;

sub _cmd {
    my $self = shift;

    return [ $self->_tool, $self->_param ];
}

sub _run_cmd {
    my $self = shift;

    $self->_open_fhs;

    # #for testing
    # say {$self->_stdout_log_fh} "running out!";
    # say {$self->_stderr_log_fh} "running err!";
    # say {$self->_cmd_log_fh} "running cmd!";
    # #for testing

    my $cmd = $self->_cmd;
    my $datestamp = DateTime->now( time_zone => 'local' ) . " ---- " . join " ", @$cmd;
    say {$self->_stdout_log_fh} $datestamp;
    say {$self->_stderr_log_fh} $datestamp;
    say {$self->_cmd_log_fh} $datestamp;
    my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $cmd, verbose => $self->verbose );
    die $error_message unless $success;
    print {$self->_stdout_log_fh} join "", @$stdout_buf;
    print {$self->_stderr_log_fh} join "", @$stderr_buf;
}

sub _open_fhs {
    my $self = shift;

    make_path( $self->_log_dir );
    open my $stdout_fh, ">>", $self->_log_dir . $self->_stdout_log;
    open my $stderr_fh, ">>", $self->_log_dir . $self->_stderr_log;
    open my $cmd_fh,    ">>", $self->_log_dir . $self->_cmd_log;
    $self->_stdout_log_fh($stdout_fh);
    $self->_stderr_log_fh($stderr_fh);
    $self->_cmd_log_fh($cmd_fh);
}

has 'out_dir' => (
    is => 'rw',
    isa => 'Str',
    traits => [qw(NoClone)],
    default => "./",
);

has '_log_dir' => (
    is      => 'rw',
    isa     => 'Str',
    traits  => [qw(NoClone)],
    lazy => 1,
    default => sub {
        my $self = shift;

        return $self->out_dir . "/logs/";
    },
);

has '_stdout_log' => (
    is => 'rw',
    isa => 'Str',
    traits => [qw(NoClone)],
    default => "/out.log",
);

has '_stderr_log' => (
    is => 'rw',
    isa => 'Str',
    traits => [qw(NoClone)],
    default => "/err.log",
);

has '_cmd_log' => (
    is => 'rw',
    isa => 'Str',
    traits => [qw(NoClone)],
    default => "/cmd.log",
);

has '_stdout_log_fh' => (
    is => 'rw',
    isa => 'FileHandle',
    traits => [qw(NoClone)],
);

has '_stderr_log_fh' => (
    is => 'rw',
    isa => 'FileHandle',
    traits => [qw(NoClone)],
);

has '_cmd_log_fh' => (
    is => 'rw',
    isa => 'FileHandle',
    traits => [qw(NoClone)],
);

has 'verbose' => (
    is      => 'rw',
    isa     => 'Bool',
);


no Moose;
__PACKAGE__->meta->make_immutable;
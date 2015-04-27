#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: sliding_read_depth.pl
#
#        USAGE: ./sliding_read_depth.pl [opts] sorted_file.bam
#
#  DESCRIPTION: Calculate depth and coverage of a sorted bamfile
#               over partitions of chroms/scaffolds with multiple output options
#
#      OPTIONS: -h --help             : show help
#               -w --windowsize <INT> : set window size
#               -m --mapped           : set windows based on mapped nucs
#               -c --compact          : do not partition into windows
#               -f --splitfile        : put chroms/scaffs in diffrent files
#               --coverage            : print mapped nucs + total nucs
#
# REQUIREMENTS: samtools, sorted bamfile as input
#         BUGS: ---
#        NOTES: Input needs to be the last argument
#       AUTHOR: Markus Lindberg  markus.lindberg89@gmail.com
# ORGANIZATION: 
#      VERSION: 1.1
#      CREATED: 19/03/2015 19:27:30 
#     REVISION: 25/03/2015 13:50:08
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;

if (scalar @ARGV < 1 ) { usage() };

foreach (@ARGV){
  if ( $_ eq "-h" || $_ eq "--help"){
    usage();
  }
}

my $FILE = pop @ARGV;

if (!-f $FILE) {
  print "$FILE is not a file\n";
  usage();
}


####=====####
#           #
#  SETUP    #
#           #
####====#####

my %opts = (
     w => 5000,
     f => 0,
);
GetOptions( \%opts,
  'w|windowsize=i',
  'm|mapped',
  'c|compact',
  'f|splitfile',
  'coverage'
) or die "Cannot parse cli arguments : $!\n";


# Hash which contains parameters and
# subrefs to test and handle data 
# from %data and %data_buffer

# Parameters
my %data_handler = (
    'handle'   => *STDOUT,         # Filehandle
    'stdout'   => 1,               # bool
    'comp'     => 'pos',           # 'pos' or 'nuc'
    'limit'    => $opts{w},        # INT
    'coverage' => $opts{coverage}, # bool
    'compact'  => $opts{c}         # bool
); 
if ($opts{m} ) { $data_handler{'comp'} = 'nuc'; };


# Methods
%data_handler = (
  %data_handler,
  'open' => sub {
    my $chr  = shift;
    my $name = shift;
    open (my $FH, '>', "${chr}_${name}.depth") or die 
        "Cannot open ${chr}_${name}.depth for writing: $!\n";

    $data_handler{'handle'} = $FH;
    $data_handler{'stdout'} = 0;
  },
  'close' => sub {
     unless ( $data_handler{'stdout'} ) {
       close $data_handler{'handle'};
       $data_handler{'handle'} = *STDOUT;
       $data_handler{'stdout'} = 1;
     }
   },
   'test' => sub {
     my $data = shift;
     my $buff = shift;

     return (($data->{$data_handler{'comp'}} - $buff->{$data_handler{'comp'}}) >= $data_handler{'limit'});
   },
   'output' => sub {
     my $data = shift;
     my $buff = shift;
     my $tot  = shift; # bool

     my $_nuc   = 0;
     my $_pos   = 0;
     my $_dep   = 0;
     my $_count = '';

     if ( $tot ) {

       $_nuc   = $data->{'nuc'};
       $_pos   = $data->{'pos'};
       $_dep   = $data->{'dep'};
       $_count = 'TOTAL';

     } else {

       $_nuc    = $data->{'nuc'} - $buff->{'nuc'};
       $_pos    = $data->{'pos'} - $buff->{'pos'};
       $_dep    = $data->{'dep'} - $buff->{'dep'};
       $_count  = $buff->{'count'};
     }
     
     my $msg = "$data->{'prev_chr'}\t$_count\t" . mean_format($_dep, $_nuc);
     if ( $data_handler{'coverage'} ) {
       $msg .= "\t$_nuc";
       $msg .= "\t$_pos";
     }

     print { $data_handler{'handle'} } "$msg\n";
   },
   'save_buffers' => sub {
     my $data = shift;
     my $buff = shift;

     map { $buff->{$_} = $data->{$_} } keys(%{$data});
     $buff->{'count'}++;
   },
   'flush' => sub {
     my $data = shift;

     $data->{'prev_chr'} = '';
     $data->{'dep'}      = 0; 
     $data->{'nuc'}      = 0;
     $data->{'pos'}      = 0;
     $data->{'count'}    = 1;
   }
);

if ($opts{c}) { $data_handler{'test'} = sub { return 0; }; }

chomp(my $st = `which samtools`);
my $prog = $st . " depth $FILE";
open my $CMD, '-|', $prog or die "Cannot open pipe $prog. $!\n exiting\n";

# Read first line to initialize data-variables
my $_firstline = <$CMD>;

my @lineparts = split /\s+/, $_firstline;

# Variables to store data
my %data = (
    'prev_chr' => $lineparts[0],
    'dep'      => $lineparts[2],
    'nuc'      => 1,
    'pos'      => $lineparts[1]
);
my %data_buffer = (
    'prev_chr' => '',
    'dep'      => 0,
    'nuc'      => 0,
    'pos'      => 0,
    'count'    => 1
);

my %total = (
    'dep'      => 0,
    'nuc'      => 0,
    'pos'      => 0
);

if ($opts{f}) {
    $data_handler{'open'}->($data{'prev_chr'}, $FILE);
}

# =============== #
#   RUN PROGRAM   #
#   ===========   #

while (my $line = <$CMD>) {
    my ($chr, $pos, $dep) = split /\s+/, $line;

    if ($data{'prev_chr'} eq $chr) {

        $data{'dep'} += $dep;
        $data{'nuc'} += 1;
        $data{'pos'}  = $pos;

        if ($data_handler{'test'}->( \%data, \%data_buffer )) {

            $data_handler{'output'}->( \%data, \%data_buffer );
            $data_handler{'save_buffers'}->( \%data, \%data_buffer );
        }

    } else {

        $data_handler{'output'}->( \%data, undef, 1 );
        $data_handler{'flush'}->( \%data_buffer );

        $total{'dep'} += $data{'dep'};
        $total{'nuc'} += $data{'nuc'};
        $total{'pos'} += $data{'pos'};

        $data{'prev_chr'} = $chr;
        $data{'dep'}      = $dep;
        $data{'nuc'}      = 1;
        $data{'pos'}      = $pos;

        if ($opts{f}) {

            $data_handler{'close'}->();
            $data_handler{'open'}->( $chr, $FILE );
        }
    }
}
close $CMD;

$data_handler{'output'}->( \%data, \%data_buffer );
$total{'prev_chr'} = 'TOTAL';
if ($opts{f}) {
    $data_handler{'close'}->();
}
$data_handler{'output'}->( \%total, undef, 1 );


#####=========#####
#
#   Subroutines
#
#==================

sub mean_format
{
  my $t = shift;
  my $n = shift;
  return(sprintf("%.3f", ($t / $n)));
}

sub usage
{
  print "
  Usage: $0 [options] <file.bam>

  Takes a _sorted_ bam-file as input and calculate read depth over the 
  chromosomes/scaffolds as they appear in the file. Need samtools in 
  \$PATH. It partitions the chromosomes/scaffolds into 'windows' and 
  calculate read depth separately on every window. The windows can either 
  be a stretch of nucleotides in the reference genome (default, suitable 
  for reseq data), or for stretches of _mapped_ nucleotides (suitable for 
  RAD data).

  Note that long stretches with unmapped nucleotides will not be partitioned
  into separate windows.

  Output (Column 4 & 5 included only with opt '--coverage'):
   Column 1: Name of chromosome/scaffold
   Column 2: Window number within chromosome/scaffold ('TOTAL' is the entire
             chromosome/scaffold)
   Column 3: Read depth
   Column 4: Number of covered nucleotides
   Column 5: Number of total nucleotides

   Options:
    -h --help               Show this help and exit
    -w --windowsize <INT>   Set the size of discrete windows (default: 5000)
    -m --mapped             Partition genome based on _mapped_ nucleotides
    -c --compact            Do not partition genome, only print stats for 
                            entire chromosomes/scaffolds
    -f --splitfile          Put results from each chromosome/scaffold into 
                            different files named <chr>_<file.bam>.depth
                            warning: many scaffolds == many files
    --coverage              Show column 4 & 5

By: Markus Lindberg, markus.lindberg89\@gmail.com\n";
  exit;
}


#!/usr/bin/env perl

## Use 'perldoc catmb.pl' or see the end of file for description.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

## Globals
my $scriptname        = $0;
my $VERSION           = '0.1.2';
my $CHANGES           = 'Mon 13 Apr 2020';
my $DEBUG             = 0;                            # 1: Debug printing
my $help              = q{};                          # q{} is empty
my $version           = q{};
my $burntrees         = q{};
my $burntrees_binary  = 'burntrees.pl';               # Change this to whatever appropriate
my $infile            = q{};
my $jump              = q{};
my $burnin            = q{};
my $pburnin           = q{};
my $burnin_arg        = q{};
my $noclose_arg       = '--noclose';
my $treesonly_arg     = '--treesonly';
my $close_arg         = '--close';
my $jump_arg          = q{};
my @output            = ();
#my $outfile = 'sdsdsd';

## Handle arguments
if (@ARGV < 1) {
    print "\n Try '$scriptname --man' for full info\n\n";
    exit(0);
}
else {
    GetOptions('help'         => sub { pod2usage(1); },
               'version'      => sub { print STDOUT "\n  $scriptname version $VERSION\n  Last changes $CHANGES\n"; exit(0) },
               'man'          => sub { pod2usage(-exitstatus => 0, -verbose => 2); },
               'jump:i'       => \$jump,
               'burnin:i'     => \$burnin,
               'pburnin:i'    => \$pburnin
              );
}

## Check if burntrees.pl can be found in path
FIND_BURNTREES:
foreach (split(/:/,$ENV{PATH})) {
    if (-x "$_/$burntrees_binary") {
        $burntrees = "$_/$burntrees_binary";
        last FIND_BURNTREES;
    }
}
if ($burntrees eq '') {
	die "Couldn't find executable $burntrees_binary (check your path or edit this script).\n\n";
}

## Create burnin argument
if ($pburnin) {
    $burnin_arg = "--pburnin=$pburnin";
}
elsif ($burnin) {
    $burnin_arg = "--burnin=$burnin";
}

## Create jump (thinning) arg
if ($jump ne '') {
    $jump_arg = "--jump=$jump";
}

## Call burntrees.pl
my $size_of_argv = @ARGV;
my $i = 0;

while (@ARGV) {
    $infile = shift(@ARGV);
    $i++;
    if ($i == 1) { # first infile
        #system("$burntrees $burnin_arg $noclose_arg $jump_arg $infile > $outfile");
        @output = `$burntrees $burnin_arg $noclose_arg $jump_arg $infile`;
    }
    elsif ($i == $size_of_argv) { # last infile
        #system("$burntrees $burnin_arg $treesonly_arg $close_arg $jump_arg $infile >> $outfile");
        my @last_output  = `$burntrees $burnin_arg $treesonly_arg $close_arg $jump_arg $infile`;
        push(@output, @last_output);
    }
    else {
        #system("$burntrees $burnin_arg $treesonly_arg $infile >> $outfile");
        my @middle_output = `$burntrees $burnin_arg $treesonly_arg $infile`;
        push(@output, @middle_output);
    }
}

## Print the output
foreach (@output) {
    print;
}


## POD documentation
=pod

=head1 NAME

catmb.pl


=head1 VERSION

Documentation for catmb.pl version 0.1.2


=head1 SYNOPSIS

catmb.pl [--burnin=<number>] [--pburnin=<number>] [--jump=<number>] FILE FILE [...] [> OUTPUT]


=head1 DESCRIPTION

Script for manipulating tree (*.t) and parameter (*.p)
files from MrBayes (v.3). This script is a helper for
concatenating several files in to one (using the same
burnin) using the script 'burntrees.pl'. burntrees.pl
needs to be installed for catmb.pl to work.


=head1 OPTIONS

Mandatory arguments to long options are mandatory for short options too


=over 8

=item B<-b, --burnin=>I<number>

Start printing after tree I<number>.


=item B<-p, --pburnin=>I<number>

Start printing after a fraction of the run, where I<number> is a percentage (e.g. "50" for half the run).


=item B<-j, --jump=>I<number>

Specify a thinning. That is, print every I<number> tree.


=item B<FILE>

Reads Nexus formatted tree B<FILEs>, preferrably a MrBayes (v.3) *.t file.
Can also read and print a MrBayes parameter file (*.p file).


=item B<OUTPUT>

Prints to B<STDOUT>.


=item B<-h, --help>

Prints help message and exits


=item B<-v, --version>

Prints version message and exits


=item B<-m, --man>

Displays the manual page


=back


=head1 USAGE

Examples:

  catmb.pl --burnin=10 run1.t run2.t run3.t > out.t
  catmb.pl --pburnin=50 run1.t run2.t run3.t > out.t
  catmb.pl --jump=10 run1.t run2.t run3.t > out.t


=head1 AUTHOR

Written by Johan A. A. Nylander


=head1 REPORTING BUGS

Please report any bugs to I<Johan.Nylander @ nrm.se>.


=head1 DEPENDENCIES

Needs burntrees.pl to run. 
Uses Perl modules Getopt::Long and Pod::Usage


=head1 DOWNLOAD

https://github.com/nylander/Burntrees


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2006-2020 Johan Nylander

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=cut


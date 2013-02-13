#!/usr/bin/perl

## Use 'perldoc burntrees.pl' or see the end of file for description.


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper;


## Globals
my $scriptname         = $0;
my $VERSION            = '0.2.0';
my $CHANGES            = '02/13/2013 10:23:38 AM';
my $DEBUG              = 0;   # Set to 1 or use --DEBUG for debug printing
my $burnin             = q{}; # q{} is empty
my $close              = q{};
#my $concatenate       = q{};
my $end                = q{};
my $figtree            = q{};
my $format             = q{};
my $getinfo            = q{};
my $help               = q{};
my $ifeellucky         = q{};
my $infile             = q{};
my $i                  = q{};
my $is_tree_file       = q{};
my $j                  = q{};
my $jump               = q{};
my $labels             = q{};
my $lastgen            = q{};
my $match              = q{};
my $nexttolastgen      = q{};
my $noclose            = q{};
my $nolabels           = q{};
my $ntrees             = q{};
my $outfile            = q{};
my $pburnin            = q{};
my $random_nr          = q{};
my $rmbrlens           = q{};
my $sci2norm           = q{};
my $start              = q{};
my $treesonly          = q{};
my $altnexus_treesonly = q{};
my $version            = q{};
my %translation_table  = ();
my $PRINT_FH;


## "MAIN" routine
MAIN:
    ## Handle arguments
    if (@ARGV < 1) {
        print "\n Try '$scriptname --man' for full info\n\n";
        exit(0);
    }
    else {
        GetOptions('help'         => sub { pod2usage(1); },
                   'version'      => sub { print STDOUT "  $scriptname version $VERSION\n  Last changes $CHANGES\n"; exit(0) },
                   'man'          => sub { pod2usage(-exitstatus => 0, -verbose => 2); },
                   'burnin:i'     => \$burnin,
                   'close'        => \$close,
                   #'concatenate'  => \$concatenate,
                   'end:i'        => \$end,
                   'format:s'     => \$format,
                   'getinfo'      => \$getinfo,
                   'ifeellucky:f' => \$ifeellucky,
                   'jump:i'       => \$jump,
                   'labels'       => \$labels,
                   'noclose'      => \$noclose,
                   'nolabels'     => \$nolabels,
                   'outfile:s'    => \$outfile,
                   'pburnin:i'    => \$pburnin,
                   'rmbrlens'     => \$rmbrlens,
                   'start:i'      => \$start,
                   'sci2norm'     => \$sci2norm,
                   'treesonly'    => \$treesonly,
                   'DEBUG'        => \$DEBUG
                  );
    }

    ## Some debug printing
    print_debug(1) if $DEBUG;

    ## Set --jump default. 1 prints every tree.
    if ($jump eq '') {
        $jump = '1';
    }

    ## --noclose is stronger than --close 
    if ($noclose) {
        $close = 0;
    }

    ## Get ONE infile from ARGV.
    if (@ARGV > 0) { 
        $infile = shift(@ARGV);
    }
    
    COUNTANDCHECK:
    ## First: count trees and check file type
    open (my $IN, '<', $infile) or die "$0 : failed to open input file $infile : $!\n";

    $ntrees = 0;
    $match = 'Peppes_Bodega';                  # Dummy match pattern to fail a match on the first lines

    while (<$IN>) {
        if (/^\s*begin\s*trees/i) {            # A .t file contains "begin trees"
            $is_tree_file = 1;
            $match = '^\s*tree';
        }
        if ((/^\s*Gen/i) or (/^\s*state/i)) {  # A .p file contains "Gen", a BEAST log file contains "STATE"
            $is_tree_file = 0;
            $match = '^\d+';
        }
        if (/$match/i) {
            $ntrees++;                         # Count trees in .t files or lines in .p files
            $nexttolastgen = $lastgen;
            if ($is_tree_file) {                          # Collect info on generation number in t file
                $_ =~ /\s+(STATE_|rep\.|gen\.)(\d+)\s+/i; # BEAST 1.7.4 | MrBayes old | MrBayes 3.2
                if ($2) {                                 # If 'rep' is not found in tree file, then
                    $lastgen = $2;
                }
                else {
                    $lastgen = 1;              # set lastgen to 1, and
                    $nexttolastgen = 0;        # nexttolastgen to 0 (to be able to get reasonable
                }                              # output from --getinfo)
            }
            else {
                $_ =~ /(^\d+)\s+/;             # Collect info on generation number in p file
                $lastgen = $1;
            }
        }
    }
    close ($IN) or warn "$0 : failed to close input file $infile : $!\n";

    ## Error if no match set (neither a .p nor a .t file)
    die "\nInfile doesn't seem to be a tree file (or a p file)\n"
        if ($match eq 'Peppes_Bodega');

    BURNIN:
    ## Set --start, --burnin, and --fburnin defaults
    if ($start eq '') {
        if ($burnin) {
            $start = $burnin + 1;
        }
        elsif ($pburnin) {
            $start =  int($pburnin * (1/100) * $ntrees) + 1;
        }
        else {
            $start = '1';
        }
    }

    ## Set probability printing defaults
    if ($ifeellucky eq '') {
        $ifeellucky = 1;
        $random_nr = 0;
    }

    ## Set --end default to $ntrees if no --end
    if ($end eq '') {
        $end = $ntrees;
    }

    BRLENS:
    ## Test brlens if rmbrlens or sci2norm
    if ($rmbrlens or $sci2norm) {
        my $has_brlens = test_brlens_presence($infile);
        if ( $has_brlens == 0 ) {
            die "\nWarning! the file doesn't seem to have branch lengths.\n\n";
        }
    }

    ## Set filehandle for printing.
    if ($outfile eq '') {
        $PRINT_FH = *STDOUT; # Using the typeglob notation in order to use STDOUT as a variable
    }
    else {
        open ($PRINT_FH, '>', $outfile) or die "$0 : Failed to open output file $outfile : $!\n\n";
    }

    FORMAT:
    ## Set output format (altnexus or phylip)
    if ($format ne '') {
        if ($is_tree_file) { # --format only applicable on tree files
            $figtree = test_figtree_format($infile); 
            if ($format =~ m/^a/i) {
                $format = 'altnexus';
                if ($nolabels) {
                    $labels = 0;
                }
                else {
                    $labels = 1;
                }
                if ($noclose) {
                    $close = 0;
                }
                if ($treesonly) {
                    $altnexus_treesonly = 1;
                }
                else {
                    print $PRINT_FH "#NEXUS\n[Trees from file $infile]\nBegin trees;\n";
                    $treesonly = 1;
                }
            }
            elsif ($format =~ m/^p/i) {
                $format = 'phylip';
                $treesonly = 1;
                $close = 0;
                if ($nolabels) {
                    $labels = 0;
                }
                else {
                    $labels = 1;
                }
            }
            else {
                die "\nWarning! output format ($format) should either be 'altnexus' or 'phylip'.\n\n"
            }
        }
        else { # Warn if --format is used with a *.p file
            die "\nWarning! argument 'format' is only applicable on tree files.\n\n"
        }
    }

    ## Read the translation table
    if ($labels) {
        %translation_table = read_labels($infile);
    }

    GETINFO:
    ## If --getinfo, print number of trees etc and quit
    if ($getinfo) {
        my $thinning = $lastgen - $nexttolastgen;
        my $totgen = $thinning * ($ntrees - 1);
        my $word = 'lines';
        if ($is_tree_file) {
            $word = 'trees';
        }
        print $PRINT_FH "Total number of $word : $ntrees.\n";
        if ($burnin or $pburnin) {
            my $Ntrees = $end - $start + 1;
            print $PRINT_FH "Number of $word after burnin : $Ntrees.\n";
        }
        if ($lastgen > 1) { # If no 'rep' in tree file, $lastgen was set to 1 (then probably not a Mrbayes *.t file)
            print $PRINT_FH "Thinning seems to have been : $thinning.\n";
            print $PRINT_FH "Number of generations was then : $totgen.\n";
        }
        if ($burnin or $pburnin) {
            my $burngen = ($start - 1) * $thinning;
            print $PRINT_FH "Burnin in generations : $burngen.\n";
        }
        if ($is_tree_file) {
            my $has_brlens = test_brlens_presence($infile);
            if ($has_brlens) {
                print $PRINT_FH "Trees appear to have branch lengths.\n";
            }
            else {
                print $PRINT_FH "Trees appear to have no branch lengths.\n";
            }
        }
        exit(0);
    }

    ## Then do some tests
    die "\nWarning! end value: $end is greater than number of trees: $ntrees.\n\n"
        if ($end > $ntrees);
    die "\nWarning! start value: $start is greater than end value: $end.\n\n"
        if ($start > $end);
    die "\nWarning! start value: $start is greater than number of trees: $ntrees.\n\n"
        if ($start > $ntrees);
    die "\nWarning! probability value: $ifeellucky is greater than one.\n\n"
        if ($ifeellucky > 1);

    ## Some debug printing
    print_debug(2) if $DEBUG;

    PRINT:
    ## Finally: do the printing
    open (my $INFILE, '<', $infile) or die "$0 : failed to open input file $infile : $!\n";
    $i = 1;
    while (<$INFILE>) {
        if ($ifeellucky < 1) { # A value less than 1 needs a random number
            $random_nr = rand;
        }
        if (/$match/i) { # If found a tree or a parameter line
            if ($sci2norm) {
                $_ = sci2norm($_); # warning, this is done one the whole string, including tree name and newick string.
            }
            #if ($rmbrlens) { # Should rmbrlens here once. No strip_brlens_print needed?
            #    #$_ = ;
            #}
            if ($format) {  # Print tree in correct format
                if ($figtree) {
                    $_ = remove_figtree_comments($_);
                }
                if ($format eq 'phylip') {
                    $_ = remove_tree_name($_);
                    if ($labels) {
                        $_ = replace_numbers($_, \%translation_table);
                    }
                }
                elsif ($format eq 'altnexus') {
                    if ($labels) {
                        $_ = replace_numbers($_, \%translation_table);
                    }
                }
                elsif ($labels) {
                    $_ = replace_numbers($_, \%translation_table);
                }
            }
            if ($i == $end) { # If last tree
                last unless $is_tree_file; ##############################################
                if ($rmbrlens) {  # If --rmbrlens,
                    print STDERR "\n=== print 1:(i:$i j:$j) ===" if $DEBUG;
                    if ($random_nr <  $ifeellucky) {
                        strip_brlens_print($_) if ( ($i % $jump) == 0 ); # Print with no branch lengths if modulus is 0.
                        if ($close) {
                            print $PRINT_FH "END;[close]\n" unless ($format eq 'phylip');
                            last;
                        }
                        elsif ($noclose) {
                            last;
                        }
                        elsif ($format eq 'altnexus') {
                            print $PRINT_FH "END;[altnex]\n" unless ($altnexus_treesonly);
                        }
                        elsif ($treesonly) {
                            last;
                        }
                        else {
                            print $PRINT_FH "END;[default]\n" unless ($format eq 'phylip');
                            last;
                        }
                    }
                }
                #elsif ($sci2norm) {
                #    ;
                #}
                else {
                    print STDERR "\n=== print 2:(i:$i j:$j) ===" if $DEBUG;
                    if ($random_nr <  $ifeellucky) {
                        print $PRINT_FH "$_" if ( ($i % $jump) == 0 ); # Print with branch lengths if modulus is 0.
                        if ($close) {
                            print $PRINT_FH "END;\n" unless ($format eq 'phylip');
                            last;
                        }
                        elsif ($noclose) {
                            last;
                        }
                        elsif ($format eq 'altnexus') {
                            print $PRINT_FH "END;\n" unless ($altnexus_treesonly);
                            last;
                        }
                        elsif ($treesonly) {
                            last;
                        }
                        else {
                            print $PRINT_FH "END;\n" unless ($format eq 'phylip');
                            last;
                        }
                    }
                }
            }
            if ($i < $end) {
                if ($i == $start) {
                    if ($rmbrlens) {
                        print STDERR "\n=== print 3:(i:$i j:$j) ===" if $DEBUG;
                        if ($random_nr <  $ifeellucky) {
                            strip_brlens_print($_); # Make sure to print the start tree.
                        }
                        $j = 0; # Set counter for --jump.
                    }
                    else {
                        print STDERR "\n=== print 4:(i:$i j:$j) ===" if $DEBUG;
                        if ($random_nr <  $ifeellucky) {
                            print $PRINT_FH "$_";
                        }
                        $j = 0;
                    }
                }
                elsif ($i > $start) {
                    if ($rmbrlens) {
                        print STDERR "\n=== print 5:(i:$i j:$j) ===" if $DEBUG;
                        if ($random_nr <  $ifeellucky) {
                            strip_brlens_print($_) if ( ($j % $jump) == 0 );
                        }
                    }
                    else {
                        print STDERR "\n=== print 6:(i:$i j:$j) ===" if $DEBUG;
                        if ($random_nr <  $ifeellucky) {
                            print $PRINT_FH "$_" if ( ($j % $jump) == 0 );
                        }
                    }
                }
            }
            $i++;
            $j++;
        }
        else { # If string is not a tree, it's either taxa descriptions or a trailing end
            if ($i > $end) {
                print STDERR "\n=== print 7:(i:$i j:$j) ===" if $DEBUG;
                print $PRINT_FH "$_" unless ($noclose);
            }
            else {
                print $PRINT_FH $_ unless ($treesonly);
            }
        }
    }
    close ($INFILE) or warn "$0 : failed to close input file $infile : $!\n";

    ## Some debug printing
    print_debug(3) if $DEBUG;

    ## Exit explicitly
    exit(0);

## End of MAIN


#===  FUNCTION  ================================================================
#         NAME:  print_debug
#      VERSION:  02/01/2013 04:30:33 PM
#  DESCRIPTION:  debug printing
#   PARAMETERS:  number
#      RETURNS:  prints to STDERR
#         TODO:  ???
#===============================================================================
sub print_debug {

    my ($number) = @_;

    print STDERR "\n\n= $number ==============\n";
    print STDERR "start:$start.\nend:$end.\nburnin:$burnin.\npburnin:$pburnin.\njump:$jump.\n";
    print STDERR "treesonly:$treesonly.\nclose:$close.\nnoclose:$noclose.\nrmbrlens:$rmbrlens.\nsci2norm:$sci2norm.\ngetinfo:$getinfo.\n";
    print STDERR "infile:$infile.\nlastgen:$lastgen.\nnexttolastgen:$nexttolastgen.\n";
    print STDERR "match:$match.\nis_tree_file:$is_tree_file.\ni:$i.\nj:$j.\nifeellucky:$ifeellucky.\nrandom_nr:$random_nr.\n";
    print STDERR "labels:$labels.\nnolabels:$nolabels.\nformat:$format.\nfigtree:$figtree.\noutfile:$outfile.\n";
    #print STDERR "concatenate: $concatenate.\n\n";
    print STDERR "==================\n";
    warn "\n (hit return to continue)\n" and getc();

} # end of print_debug


#===  FUNCTION  ================================================================
#         NAME:  read_labels
#      VERSION:  02/03/2013 01:07:37 PM
#  DESCRIPTION:  reads sequence (taxon) labels and associates them with numbers
#   PARAMETERS:  file name
#      RETURNS:  Hash with taxon translation table
#         TODO:  ???
#===============================================================================
sub read_labels {

    my ($file) = @_;
    my %hash = (); # key: number, value: name

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        my $line = $_;
        chomp($line);
        if($line =~ m/\s*tree\s+state/i) {
            last;
        }
        elsif($line =~ m/^\s*(\d+)\s+([\w|\s|\W]+)$/) { # capture the number, and the
            my $number = $1;                            # taxon name allowing for single-
            my $name = $2;                              # quoted taxon names
            $name =~ s/\s*,\s*$//;
            $name =~ s/\s*;\s*$//;
            $hash{$number} = $name;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    if ($DEBUG) {
        print STDERR "\nhash in read_labels:\n";
        for my $key ( keys %hash ) {
            my $value = $hash{$key};
            print STDERR "$key => $value\n";
        }
    }

    return %hash;

} # end of read_labels

## Keep this alternative code for reading labels from MrBayes/BEAST?
#sub read_labels {
#
#    my ($file) = @_;
#    my %hash = ();
#
#    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
#    while(<$FILE>) {
#        if((/\s*tree\s+rep/i) or ((/\s*tree\s+state/i))) {
#            last;
#        }
#        if(/^\s+\d+\s+\b/) {
#            my ($number, $name_x) = split;
#            my @letters = split //, $name_x;
#            if ( ($letters[-1] eq ",") or ($letters[-1] eq ";") ) {
#                pop @letters;
#            }
#            my $name = join "", @letters;
#            $hash{$number} = $name;
#        }
#    }
#    close ($FILE) or warn "$0 : failed to close file $file : $!\n";
#
#    if ($DEBUG) {
#        print STDERR "\nhash in read_labels:\n";
#        for my $key ( keys %hash ) {
#            my $value = $hash{$key};
#            print STDERR "$key => $value\n";
#        }
#    }
#
#    return %hash;
#
#} # end of read_labels


#===  FUNCTION  ================================================================
#         NAME:  remove_figtree_comments
#      VERSION:  01/25/2013 09:34:53 PM
#  DESCRIPTION:  Removes the figtree comments
#   PARAMETERS:  tree string
#      RETURNS:  tree string without figtree comments
#         TODO:  ?
#===============================================================================
sub remove_figtree_comments {

    my ($figtree) = @_;

    $figtree =~ s/(\[.+?\])//g;

    return($figtree);

} # end of remove_figtree_comments


#===  FUNCTION  ================================================================
#         NAME:  remove_tree_name
#      VERSION:  02/03/2007 02:31:54 AM CET
#  DESCRIPTION:  Removes the tree name
#   PARAMETERS:  Line with tree string. Uses global variable $infile
#      RETURNS:  Tree string
#         TODO:  The regexp will not work correctly if there is a space between
#                the last closing bracket and the trailing ';'. Could be solved
#                like this:
#                if ( $piece =~ ^\(.+\)$ ) { # starts with ''( and ends in ')'
#                    $tree = $piece . ";";
#                }
#                elsif ($piece =~ ^\(.+[\);]$) { # starts with ''( and ends in ');'
#                    $tree = $piece . ";";
#                }
#                Any spaces in the tree description will of course create
#                truncated trees, however.
#===============================================================================
sub remove_tree_name {
    
    my ($line) = @_;
    my $tree   = '';

    my @pieces = split /\s+/, $line;

    foreach my $piece (@pieces) {
        if ($piece =~ /^\(.+[\);]$/) { # if piece starts with '(' and ends in ');'
            $tree = "$piece\n";        # hack: needed to add line break
        }
    }

    if ($tree eq '') {
        die "Warning: Could not read tree format.\nAcceptable format: tree name = (1,(2,3));\n";
    }

    return $tree;

} # end of remove_tree_name


#===  FUNCTION  ================================================================
#         NAME:  replace_numbers
#      VERSION:  02/03/2013 04:13:04 PM
#  DESCRIPTION:  replaces numbers with sequence (taxon) labels.
#                Possible matches and replacements are
#                    ',123:' => ',foo:'
#                    ',123,' => ',foo,'
#                    ',123)' => ',foo)'
#                    ',123[' => ',foo['
#                    '(123:' => '(foo:'
#                    '(123,' => '(foo,'
#                    '(123[' => '(foo['
#   PARAMETERS:  tree string and reference to hash holding translation table.
#                Uses global variable $infile
#      RETURNS:  tree string with seq labels instead of numbers
#         TODO:  Test FigTree format
#===============================================================================
sub replace_numbers {

    my ($tree, $hash_reference) = @_;

    foreach my $number (keys %$hash_reference) { # $key is number, $value is label

        my $label = $hash_reference->{$number};

        if ( $tree =~ /([,|\(])$number([\)|,|:|\[])/ ) {
            $tree =~ s/([,|\(])$number([\)|,|:|\[])/$1$label$2/;
        }
    }

    return $tree;

} # end of replace_numbers


#===  FUNCTION  ================================================================
#         NAME: sci2norm
#      VERSION: 02/03/2013 05:44:21 PM
#  DESCRIPTION: translate tree string with branch lengths in scientific notation
#               to normal notation
#   PARAMETERS: tree string
#      RETURNS: tree string
#         TODO: If FigTree format, translate all scientific notations, not only brlens.
#               First, have this function translate all scientific notations, not only branch lengths.
#===============================================================================
sub sci2norm {

    my ($tree) = @_;

    #$tree =~ s/:(\d+\.\d+e[\+|\-]\d+)/sci2norm_colon($1)/ieg;
    $tree =~ s/(\d+\.\d+e[\+|\-]\d+)/sci2norm_no_colon($1)/ieg;

    return($tree);

} # end of sci2norm 


#===  FUNCTION  ================================================================
#         NAME: sci2norm_print
#      VERSION: 02/08/2013 03:16:45 PM
#  DESCRIPTION: translate string of scientific number notation to decimal number
#   PARAMETERS: ????
#      RETURNS: string (decimal number) "89.001"
#         TODO: ????
#===============================================================================
sub sci2norm_print {

    my($sci) = @_;

    my $nr = sprintf "%f", $sci;

    return($nr);

} # end of sci2norm_print


#===  FUNCTION  ================================================================
#         NAME:  strip_brlens_print
#      VERSION:  02/13/2013 09:25:00 AM
#  DESCRIPTION:  Removes the branch lengths from a tree descriptions and print
#   PARAMETERS:  tree string
#      RETURNS:  Void. Prints to $PRINT_FH
#         TODO:  ?
#===============================================================================
sub strip_brlens_print {

    ($_) = @_; 
    
    if (/e-\d+/i) { # if scientific notation ":1.309506485347851e-01" or ":1.309506485347851E-01"
        $_ =~ s/:[\d\.e\-]+//gi;
    }
    else {
        $_ =~ s/:[\d\.]+//g; # Search for any pattern such as ":0.033372" and replace with nothing
    }

    print $PRINT_FH $_;

} # end of strip_brlens_print


#===  FUNCTION  ================================================================
#         NAME:  test_brlens_presence
#      VERSION:  02/02/2007 12:12:05 AM CET
#  DESCRIPTION:  Tests for presence of branch lengths in the tree description.
#                Warning: no error checking. Assumes Nexus tree format.
#   PARAMETERS:  string containing tree description
#      RETURNS:  1: brlens present, 0: brlens absent
#         TODO:  ?
#===============================================================================
sub test_brlens_presence {

    my ($file) = @_;
    my $brl    = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        if(/^\s*tree/i) {
            if ( $_ =~ /:/ ) {
                $brl = 1;
            }
            last;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    return $brl;

} # end of test_brlens_presence


#===  FUNCTION  ================================================================
#         NAME:  test_figtree_format
#      VERSION:  02/02/2007 12:12:05 AM CET
#  DESCRIPTION:  Tests for figtree format, i.e., presence comments in square
#                brackets within the tree string
#                Warning: no extensive error checking. Looks only for 'tree' and '&'
#   PARAMETERS:  string containing tree description
#      RETURNS:  1: figtree format, 0: not figtree format
#         TODO:  ?
#===============================================================================
sub test_figtree_format {

    my ($file)  = @_;
    my $figtree = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        if(/^\s*tree/i) {
            if ( $_ =~ /&/ ) { # Presence of '&' probably is a figtree tree
                $figtree = 1;
            }
            last;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    return $figtree;

} # end of test_figtree_format




#===  POD DOCUMENTATION  =======================================================
#      VERSION:  02/13/2013 10:17:54 AM
#  DESCRIPTION:  Documentation
#         TODO:  ?
#===============================================================================
POD:
=pod

=head1 NAME

burntrees.pl


=head1 VERSION

Documentation for burntrees.pl version 0.2.0


=head1 SYNOPSIS

burntrees.pl [--burnin=<number>] [--pburnin=<number>] [--start=<number>] [--end=<number>] [--jump=<number>] [--IFeelLucky=<number>] [--treesonly] [--rmbrlens] [--sci2norm] [--[no]close] [--getinfo] [--[no]labels] [--format=altnexus|phylip] [--outfile=<file_name>] FILE [> OUTPUT]


=head1 DESCRIPTION

Script for manipulating tree (*.t, *.trprobs, *.con, *.trees) and parameter (*.p, *.log) files
from MrBayes (v.3) or BEAST.

The script extracts trees and (by default) the taxon translation table and the trailing "end;"
from tree file.

A number of options are available:

Any contiguous interval of trees can be printed, as well as trees only (nothing other than tree descriptions).

The samples can be thinned by setting a value for how many trees to jump before next is printed.

Branch lengths (if present) can be removed from trees before printing.

Branch lengths in scientific numeric format can be transformed to a fixed numeric format.

A random set of trees can be printed from the tree file.

Lines can also be extracted from a MrBayes *.p file.

Trees can be printed in Phylip (Newick) format or as altnexus (sequence labels instead of numbers),
that is, the script can serve as an efficient tree format converter.



=head1 OPTIONS

Mandatory arguments to long options are mandatory for short options too


=over 8

=item B<-b, --burnin=>I<number>

Start printing after tree I<number>.


=item B<-c, --close>

Forces a trailing "end;" to be printed after the last tree.
B<--noclose> prevents the "end;" to be printed.
Note that the trailing "end;" in the tree file is printed
by default unless B<--noclose> is given.


=item B<-co, --concatenate>

Concatenate several files -- I<Not yet implemented>. See B<USAGE> for alternatives.


=item B<-e, --end=>I<number>

End the printing of trees after tree I<number> (inclusively).
If no B<--end> is given, prints to last tree in file.


=item B<-f, --format=>I<format>

Trees are printed as specified by I<format>, where I<format> is either I<altnexus>: with sequence (taxon) labels instead of numbers, or I<phylip> (the Newick format). 


=item B<-g, --getinfo>

Print information about the number of trees (or samples in p file), thinning and number of samples in file and quit.


=item B<-h, --help>

Prints help message and exits.


=item B<-i, --IFeelLucky=>I<number>

Specify a probability (value between 0 -- 1) for each tree to be printed. That is, print each tree with prob. I<number>.
Note that B<--IFeelLucky> has precedence over B<--jump>.


=item B<-j, --jump=>I<number>

Specify a thinning. That is, print every I<number> tree.


=item B<-l, --labels>

Print trees using sequence (taxon) labels instead of the sequence numbers from the translation table.
B<--nolabels> (which is the default) prevents the sequence numbers to be substituted.


=item B<-m, --man>

Displays the manual page.


=item B<-o, --outfile=>I<file_name>

Print directly to file I<file_name> instead of standard out.


=item B<-p, --pburnin=>I<number>

Start printing after a fraction of the run, where I<number> is a percentage (e.g. "50" for half the run).


=item B<-r, --rmbrlens>

Remove branch lengths from trees.


=item B<-sc, --sci2norm>

Translate branch lengths from scientific to normal or fixed.


=item B<-st, --start=>I<number>

Start printing from tree I<number> (inclusively).


=item B<-t, --treesonly>

Print trees only (do not print taxon descriptions etc.).
If used on a *.p file it skips the "ID" line and the headers.


=item B<-v, --version>

Prints version message and exits.


=item B<FILE>

Reads a Nexus formatted tree B<FILE>, preferrably MrBayes (v.3) *.t,
*.con, and *.trprobs files. B<FILE> can also be a MrBayes parameter
file (*.p), or any Nexus formatted tree file, e.g., output from
BEAST or other phylogenetic MCMC software (not thoroughly tested!).


=item B<OUTPUT>

Prints to B<STDOUT> unless B<--outfile=> is used.


=back


=head1 USAGE

Examples:

  burntrees.pl --burnin=10 data.t > out.t
  burntrees.pl --pburnin=50 data.t
  burntrees.pl --start=11 --end=30 data.t
  burntrees.pl --jump=10 data.t
  burntrees.pl --treesonly data.t
  burntrees.pl --getinfo -b=10 data.t
  burntrees.pl --rmbrlens data.t
  burntrees.pl --ifeellucky=0.50 data.t
  burntrees.pl -b=10 -j=10 -t -r data.t
  burntrees.pl --treesonly -b=10 data.p
  burntrees.pl --format=altnexus data.t
  burntrees.pl -f=phylip --outfile=data.phy data.t
  burntrees.pl -f=p -b=1 data.con


For concatenation of several files, use (note the append redirection, ">>". See also separate script catmb.pl):

  burntrees.pl -b=10 -noc  data.run1.t >  data.t
  burntrees.pl -b=40 -t    data.run2.t >> data.t
  burntrees.pl -b=20 -t -c data.run3.t >> data.t


To print the '#NEXUS', 'begin trees;', and the translation table only (no trees), use

  burntrees.pl -i=0 -noc data.t


To convert the MrBayes file to a 'altnexus' file or a 'phylip' file, use

  burntrees.pl --format=altnexus data.t
  burntrees.pl --format=phylip data.t


To concatenate several files in to one altnexus formatted file, use (note the combination of '--format' and '--labels')

  burntrees.pl -b=10 -noc  -f=a data.run1.t >  data.t
  burntrees.pl -b=40 -t    -l   data.run2.t >> data.t
  burntrees.pl -b=20 -t -c -l   data.run3.t >> data.t


To extract the second tree in the MrBayes *.con file in phylip format, use

  burntrees.pl -b=1 -f=p data.con


To change the branch length format from scientific to numerical use

  burntrees.pl --sci2norm data.con.tre

 

=head1 AUTHOR

Written by Johan A. A. Nylander


=head1 REPORTING BUGS

Please report any bugs to I<Johan.Nylander @ bils.se>.


=head1 DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage


=head1 DOWNLOAD

https://github.com/nylander/Burntrees


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2006--2013, Johan Nylander.
All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 
http://www.gnu.org/copyleft/gpl.html 


=cut


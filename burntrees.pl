#!/usr/bin/env perl

## Use 'perldoc burntrees.pl' or see the end of file for description.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


## Globals
my $scriptname         = $0;
my $VERSION            = '0.3.1';
my $CHANGES            = 'Mon 13 Apr 2020';
my $DEBUG              = 0;   # Set to 1 or use --DEBUG for debug printing
my $burnin             = q{};
my $close              = q{};
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
my $myr                = q{};
my $nexttolastgen      = q{};
my $noclose            = q{};
my $nolabels           = q{};
my $ntrees             = q{};
my $outfile            = q{};
my $pburnin            = q{};
my $random_nr          = q{};
my $rmbrlens           = q{};
my $rmcomments         = q{};
my $rmsupport          = q{};
my $sci2norm           = q{};
my $seed               = q{};
my $start              = q{};
my $treesonly          = q{};
my $altnexus_treesonly = q{};
my $version            = q{};
my $was_newick         = q{};
my %translation_table  = ();
my $PRINT_FH;


## "MAIN" routine
MAIN:
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
                   'DEBUG+'       => \$DEBUG,
                   'end:i'        => \$end,
                   'format:s'     => \$format,
                   'getinfo'      => \$getinfo,
                   'ifeellucky:f' => \$ifeellucky,
                   'jump:i'       => \$jump,
                   'labels'       => \$labels,
                   'myr'          => \$myr,
                   'noclose'      => \$noclose,
                   'nolabels'     => \$nolabels,
                   'outfile:s'    => \$outfile,
                   'pburnin:i'    => \$pburnin,
                   'rmbrlens'     => \$rmbrlens,
                   'rmcomments'   => \$rmcomments,
                   'rmsupport'    => \$rmsupport,
                   'sci2norm:-1'  => \$sci2norm,
                   'seed:i'       => \$seed,
                   'start:i'      => \$start,
                   'treesonly'    => \$treesonly,
                  );
    }

    ## Some debug printing
    print_debug(1) if ($DEBUG > 1);

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
        chomp;
        if (/$match/i) {
            $ntrees++;                         # Count trees in .t files or lines in .p files
            $nexttolastgen = $lastgen;
            if ($is_tree_file) {                               # Collect info on generation number in t file
                $_ =~ /\s+(STATE_|rep\.|gen\.)(\d+)[\s+|\[]/i; # BEAST 1.8.3 | MrBayes old | MrBayes 3.2
                if ($2) {                                      # If 'rep' is not found in tree file, then
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
        elsif (/^\s*begin\s+trees/i) {         # A .t file contains "begin trees"
            $is_tree_file = 1;
            $match = '^\s*tree';
            next;
        }
        elsif (/\)\s*;\s*$/) {                 # A phylobase tree file have newick trees on one line
            $is_tree_file = 1;
            $was_newick = 1;
            $match = '^\s*\(';
            $ntrees++;
            next;
        }
        elsif ((/^\s*Gen/i) or (/^\s*state/i) or (/^\s*#cycle/i) ) { # A .p file contains "Gen", a BEAST log file contains "STATE",
            $is_tree_file = 0;                                       # A phylobase .trace file contains "#cycle"
            $match = '^\s*\d+';
            next;
        }
    }
    close ($IN) or warn "$0 : failed to close input file $infile : $!\n";

    ## Error if no match set (neither a .p nor a .t file)
    if ($match eq 'Peppes_Bodega') {
        die "\nInfile doesn't seem to be a tree file (or a p file)?\n"
    };

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
    elsif ($seed) {
        srand $seed;
    }

    ## Set --end default to $ntrees if no --end
    if ($end eq '') {
        $end = $ntrees;
    }

    BRLENS:
    ## Test brlens if rmbrlens
    if ($rmbrlens) {
        my $has_brlens = test_has_brlens($infile);
        if ( $has_brlens == 0 ) {
            die "\nWarning! the file doesn't seem to have branch lengths.\n\n";
        }
    }

    ## Test if support values
    if ($rmsupport) {
        my $has_support = test_has_support($infile);
        if ( $has_support == 0 ) {
            die "\nWarning! the file doesn't seem to have support values.\n\n";
        }
    }

    ## Test clockrate if myr
    if ($myr) {
        my $has_clockrate = test_has_clockrate($infile);
        if ( $has_clockrate == 0 ) {
            die "\nWarning! the file doesn't seem to have clock rates.\n\n";
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
        my $thinning;
        my $totgen;
        if ($lastgen) {
            $thinning = $lastgen - $nexttolastgen;
            $totgen = $thinning * ($ntrees - 1);
        }
        my $word = 'lines';
        if ($is_tree_file) {
            $word = 'trees';
        }
        print $PRINT_FH "Total number of $word : $ntrees\n";
        if ($burnin or $pburnin) {
            my $Ntrees = $end - $start + 1;
            print $PRINT_FH "Number of $word after burnin : $Ntrees\n";
        }
        if (($lastgen) && ($lastgen > 1)) { # If no 'rep' in tree file, $lastgen was set to 1 (then probably not a Mrbayes *.t file)
            print $PRINT_FH "Thinning seems to have been : $thinning\n";
            print $PRINT_FH "Number of generations was then : $totgen\n";
        }
        if ($burnin or $pburnin) {
            my $burngen = ($start - 1) * $thinning;
            print $PRINT_FH "Burnin in generations : $burngen\n";
        }
        if ($is_tree_file) {
            my $has_brlens = test_has_brlens($infile);
            if ($has_brlens) {
                print $PRINT_FH "Trees appear to have branch lengths.\n";
            }
            else {
                print $PRINT_FH "Trees appear to have no branch lengths.\n";
            }
            my $has_support = test_has_support($infile);
            if ($has_support) {
                print $PRINT_FH "Trees appear to have support values.\n";
            }
            else {
                print $PRINT_FH "Trees appear to have no support values.\n";
            }
            my $is_figtree = test_figtree_format($infile);
            if ($is_figtree) {
                print $PRINT_FH "Trees appear to be in figtree format.\n";
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
    print_debug(2) if ($DEBUG > 1);

    PRINT:
    ## Finally: do the printing
    open (my $INFILE, '<', $infile) or die "$0 : failed to open input file $infile : $!\n";
    $i = 1;
    while (<$INFILE>) {
        chomp;
        if (/$match/i) { # If found a tree or a parameter line
            if ($ifeellucky < 1) { # A value less than 1 needs a random number
                $random_nr = rand;
            }
            if ($sci2norm) {
                $_ = sci2norm($_,$sci2norm); # warning, this is done one the whole string, including tree name and newick string.
            }
            if ($myr) {
                $_ = brlen2time($_);
            }
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
                    if ($was_newick) {
                        $_ = "  tree rep.$i = " . $_;
                    }
                }
                elsif ($labels) {
                    $_ = replace_numbers($_, \%translation_table);
                }
            }
            if ($rmcomments) {
                $_ = remove_figtree_comments($_);
            }
            if ($rmsupport) {
                $_ = remove_support($_);
            }
            if ($i == $end) { # If last tree
                if ($rmbrlens) {
                    print STDERR "\n=== print 1:(i:$i j:$j) ===" if $DEBUG;
                    if ($random_nr <  $ifeellucky) {
                        strip_brlens_print($_) if ( ($i % $jump) == 0 ); # Print with no branch lengths if modulus is 0.
                        if ($close) {
                            my $l = $DEBUG ? "END; [DEBUG:close]\n" : "END;\n";
                            print $PRINT_FH $l unless ($format eq 'phylip');
                            last;
                        }
                        elsif ($noclose) {
                            last;
                        }
                        elsif ($format eq 'altnexus') {
                            my $l = $DEBUG ? "END; [DEBUG:altnex]\n" : "END;\n";
                            print $PRINT_FH $l unless ($altnexus_treesonly);
                            last;
                        }
                        elsif ($treesonly) {
                            last;
                        }
                        else {
                            my $l = $DEBUG ? "END; [DEBUG:default]\n" : "END;\n";
                            print $PRINT_FH $l unless ($format eq 'phylip');
                            last;
                        }
                    }
                }
                else {
                    print STDERR "\n=== print 2:(i:$i j:$j) ===" if $DEBUG;
                    if ($random_nr <  $ifeellucky) {
                        my $l = $DEBUG ? "$_ [DEBUG:bpa]\n" : "$_\n";
                        print $PRINT_FH $l if ( ($i % $jump) == 0 ); # Print with branch lengths if modulus is 0.
                        if ($close) {
                            my $l = $DEBUG ? "END; [DEBUG:cpa]\n" : "END;\n";
                            print $PRINT_FH $l unless ($format eq 'phylip');
                            last;
                        }
                        elsif ($noclose) {
                            last;
                        }
                        elsif ($format eq 'altnexus') {
                            my $l = $DEBUG ? "END; [DEBUG:dpa]\n" : "END;\n";
                            print $PRINT_FH $l unless ($altnexus_treesonly);
                            last;
                        }
                        elsif ($treesonly) {
                            last;
                        }
                        elsif ($is_tree_file) {
                            my $l = $DEBUG ? "END; [DEBUG:epa]\n" : "END;\n";
                            print $PRINT_FH $l unless ($format eq 'phylip');
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
                        $j = 1;
                    }
                    else {
                        print STDERR "\n=== print 4:(i:$i j:$j) ===" if $DEBUG;
                        if ($random_nr <  $ifeellucky) {
                            my $l = $DEBUG ? "$_ [DEBUG:fpa]\n" : "$_\n";
                            print $PRINT_FH $l;
                        }
                        $j = 1;
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
                            my $l = $DEBUG ? "$_ [DEBUG:gpa]\n" : "$_\n";
                            print $PRINT_FH $l if ( ($j % $jump) == 0 );
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
                my $l = $DEBUG ? "$_ [DEBUG:hpa]\n" : "$_\n";
                print $PRINT_FH $l unless ($noclose);
            }
            else {
                my $l = $DEBUG ? "$_ [DEBUG:ipa]\n" : "$_\n";
                print $PRINT_FH $l unless ($treesonly);
            }
        }
    }
    close ($INFILE) or warn "$0 : failed to close input file $infile : $!\n";

    ## Some debug printing
    print_debug(3) if ($DEBUG > 1);

    ## Exit explicitly
    exit(0);

## End of MAIN


#===  FUNCTION  ================================================================
#         NAME: brlen2time
#      VERSION: 04/09/2013 01:49:02 PM
#  DESCRIPTION: Transform branch lengths from substitutions per site to time units
#   PARAMETERS: tree string
#   tree gen.5000000[&B Igrbranchlens{all}] = [&R] [&clockrate=8.332129945298162e-04] (1:4.907239999154218e-02[&B Igrbranchlens{all} 5.967617061002457e-02],(...
#      RETURNS:  tree string
#===============================================================================
sub brlen2time {

    my ($tree) = @_;

    $tree =~ m/&clockrate\s*=\s*(\d+\.?\d*([eE][+-]?\d+)?)\]/i;
    my $clockrate = $1;

    ## Iterate over tree and substitute all branch lengths with the age:
    $tree =~ s/:(\d+\.?\d*([eE][+-]?\d+)?)/brlen2time_print($clockrate,$1)/ieg;

    return($tree);

} # end brlen2time


#===  FUNCTION  ================================================================
#         NAME: brlen2time_print
#      VERSION: 04/10/2013 09:34:53 AM
#  DESCRIPTION: print branch lengths in time units PLUS an ad hoc ':'!
#   PARAMETERS: clockrate, branch length
#      RETURNS: string (decimal number) and a colon(!) :89.001
#         TODO: Work around the ':' when used with brlen2time
#===============================================================================
sub brlen2time_print {

    my($cr, $bl) = @_;

    my $myr = sprintf ":%f", $bl/$cr;

    return($myr);

} # end of brlen2time_print


#===  FUNCTION  ================================================================
#         NAME: fig2simple_with_pp
#      VERSION: 08/11/2016 03:10:48 PM
#  DESCRIPTION: Transform figtree format to "simple" format with posterior probabilities
#   PARAMETERS: tree in figtree format with scientific dec format
#      RETURNS: newick string in "simple" format
#         TODO: Check if sci format?
#               Testing!
#               Note: look at the possibility to have the "--format" argument handle
#               the figtree to simple conversion!!!
#===============================================================================
#sub fig2simple_with_pp {
#
#    my ($tree) = @_;
#
#    if ($tree =~ /&prob=/i) {
#        $tree =~ s/\)\[&prob=(\d+\.\d+e[+-]?\d+)[^\]]+\]?/\)$1/g;
#    }
#
#    $tree = remove_figtree_comments($tree);
#
#    #$tree = sci2norm($tree);
#
#    return($tree);
#
#} # end of fig2simple_with_pp


#===  FUNCTION  ================================================================
#         NAME: print_debug
#      VERSION: 08/22/2016 02:01:53 PM
#  DESCRIPTION: debug printing
#   PARAMETERS: number
#      RETURNS: prints to STDERR
#===============================================================================
sub print_debug {

    my ($number) = @_;

    print STDERR "\n\n= $number ==============\n";
    print STDERR "burnin:$burnin.\n";
    print STDERR "close:$close.\n";
    print STDERR "end:$end.\n";
    print STDERR "figtree:$figtree.\n";
    print STDERR "format:$format.\n";
    print STDERR "getinfo:$getinfo.\n";
    print STDERR "help:$help.\n";
    print STDERR "i:$i.\n";
    print STDERR "ifeellucky:$ifeellucky.\n";
    print STDERR "infile:$infile.\n";
    print STDERR "is_tree_file:$is_tree_file.\n";
    print STDERR "j:$j.\n";
    print STDERR "jump:$jump.\n";
    print STDERR "labels:$labels.\n";
    print STDERR "lastgen:$lastgen.\n";
    print STDERR "match:$match.\n";
    print STDERR "myr:$myr.\n";
    print STDERR "nexttolastgen:$nexttolastgen.\n";
    print STDERR "noclose:$noclose.\n";
    print STDERR "nolabels:$nolabels.\n";
    print STDERR "ntrees:$ntrees.\n";
    print STDERR "outfile:$outfile.\n";
    print STDERR "pburnin:$pburnin.\n";
    print STDERR "random_nr:$random_nr.\n";
    print STDERR "rmbrlens:$rmbrlens.\n";
    print STDERR "rmcomments:$rmcomments.\n";
    print STDERR "rmsupport:$rmsupport.\n";
    print STDERR "sci2norm:$sci2norm.\n";
    print STDERR "seed:$seed.\n";
    print STDERR "start:$start.\n";
    print STDERR "treesonly:$treesonly.\n";
    print STDERR "altnexus_treesonly:$altnexus_treesonly.\n";
    print STDERR "==================\n";

    warn "\n (hit return to continue)\n" and getc();

} # end of print_debug


#===  FUNCTION  ================================================================
#         NAME: read_labels
#      VERSION: 02/03/2013 01:07:37 PM
#  DESCRIPTION: reads sequence (taxon) labels and associates them with numbers
#   PARAMETERS: file name
#      RETURNS: Hash with taxon translation table
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


#===  FUNCTION  ================================================================
#         NAME: remove_figtree_comments
#      VERSION: 08/12/2016 01:15:23 PM
#  DESCRIPTION: Removes the figtree comments. If "&prob" are present in the
#               description (as in MrBayes con.tre files), the probabilities
#               are extracted and included in the newick string!
#               Currently (mb3.2.7), prob numbers are in scientific format.
#   PARAMETERS: tree string
#      RETURNS: tree string without figtree comments
#===============================================================================
sub remove_figtree_comments {

    my ($tree) = @_;

    if ($tree =~ /&prob=\d+\.\d+e[+-]?/i) {
        $tree =~ s/\)\[&prob=(\d+\.\d+e[+-]?\d+)[^\]]+\]?/\)$1/ig;
    }
    elsif ($tree =~ /&prob=\d+\.\d+[^eE]?/i) {
        $tree =~ s/\)\[&prob=(\d+\.\d+?)[^\]]+\]?/\)$1/ig;
    }

    $tree =~ s/(\[.+?\])//g;

    return($tree);

} # end of remove_figtree_comments


#===  FUNCTION  ================================================================
#         NAME: remove_support
#      VERSION: 08/16/2016 01:20:53 PM
#  DESCRIPTION: Removes the support values from simple newick string:
#               )1.00000000e+00) -> ))
#               )1.00000000e+00:3.231983e-02) -> ):3.231983e-02)
#               )100) -> ))
#               )100:0.001) -> ):0.001)
#
#   PARAMETERS: tree string
#      RETURNS: tree string without support values
#         TODO: testing!
#===============================================================================
sub remove_support {

    my ($tree) = @_;

    $tree =~ s/\)[\d\.Ee+-]+/\)/g;

    return($tree);

} # end of remove_support


#===  FUNCTION  ================================================================
#         NAME: remove_tree_name
#      VERSION: 08/29/2016 11:20:18 AM
#  DESCRIPTION: Removes the tree name
#   PARAMETERS: Line with tree string. Uses global variable $infile
#      RETURNS: Tree string
#         TODO: Any spaces in the tree description will create
#               truncated trees!
#               Make sure the hack is no longer needed!
#===============================================================================
sub remove_tree_name {

    my ($line) = (@_);
    my $tree   = '';

    my @pieces = split /\s+/, $line;

    foreach my $piece (@pieces) {
        if ($piece =~ /^\([^;]+[\);]$/) { # if piece starts with '(' and ends in ')' or ';'
            #$tree = "$piece\n";       # hack: needed to add line break
            $tree = "$piece";
        }
    }

    if ($tree eq '') {
        die "Warning: Could not read tree format.\nAcceptable format: tree name = (1,(2,3));\n";
    }

    return $tree;

} # end of remove_tree_name


#===  FUNCTION  ================================================================
#         NAME: replace_numbers
#      VERSION: 02/03/2013 04:13:04 PM
#  DESCRIPTION: Replaces numbers with sequence (taxon) labels.
#               Possible matches and replacements are
#                   ',123:' => ',foo:'
#                   ',123,' => ',foo,'
#                   ',123)' => ',foo)'
#                   ',123[' => ',foo['
#                   '(123:' => '(foo:'
#                   '(123,' => '(foo,'
#                   '(123[' => '(foo['
#   PARAMETERS: tree string and reference to hash holding translation table.
#               Uses global variable $infile
#      RETURNS: tree string with seq labels instead of numbers
#         TODO: Test FigTree format
#===============================================================================
sub replace_numbers {

    my ($tree, $hash_reference) = @_;

    foreach my $number (keys %$hash_reference) { # $key is number, $value is label

        my $label = $hash_reference->{$number};

        if ( $tree =~ /([,\(])$number([\),:\[])/ ) {
            $tree =~ s/([,\(])$number([\),:\[])/$1$label$2/;
        }
    }

    return $tree;

} # end of replace_numbers


#===  FUNCTION  ================================================================
#         NAME: sci2norm
#      VERSION: 07/10/2014 04:16:09 PM
#  DESCRIPTION: Translate tree string with branch lengths in scientific notation
#               to normal notation.
#   PARAMETERS: tree string
#      RETURNS: tree string
#         TODO: test on p files with scientific notation
#===============================================================================
sub sci2norm {

    my ($tree,$dec) = (@_);

    $tree =~ s/(\d+\.\d+e[+-]?\d+)/sci2norm_print($1,$dec)/ieg;

    return($tree);

} # end of sci2norm


#===  FUNCTION  ================================================================
#         NAME: sci2norm_print
#      VERSION: 08/12/2016 05:36:38 PM
#  DESCRIPTION: Translate string of scientific number notation to decimal number.
#   PARAMETERS: string (scientific notation, e.g. "1.600103158188143e-02")
#      RETURNS: string (decimal number, e.g. "0.016001")
#===============================================================================
sub sci2norm_print {

    my($sci,$dec) = (@_);
    my $nr;
    if ($dec == -1) {
        $nr = sprintf "%f", $sci;
    }
    else {
        $nr = sprintf "%.${dec}f", $sci;
    }

    return($nr);

} # end of sci2norm_print


#===  FUNCTION  ================================================================
#         NAME: strip_brlens_print
#      VERSION: Mon 13 Apr 2020 12:59:44 PM CEST
#  DESCRIPTION: Removes the branch lengths from a tree descriptions and print.
#   PARAMETERS: tree string
#      RETURNS: Void. Prints to $PRINT_FH
#===============================================================================
sub strip_brlens_print {

    ($_) = @_;

    if (/e[-+]\d+/i) { # ":1.309506485347851e-01" or ":1.309506485347851E-01" or ":1.053210e+00"
        $_ =~ s/:[e\d\.\-\+]+//gi;
    }
    else {
        $_ =~ s/:[\d\.]+//g; # Search for any pattern such as ":0.033372" and replace with nothing
    }

    print $PRINT_FH "$_\n";

} # end of strip_brlens_print


#===  FUNCTION  ================================================================
#         NAME: test_has_brlens
#      VERSION: 09/29/2016 05:31:03 PM
#  DESCRIPTION: Tests for presence of branch lengths in the tree description.
#               Warning: no error checking. Assumes Nexus or Newick tree format,
#               and tree description without linebreaks.
#   PARAMETERS: string containing tree description
#      RETURNS: 1: brlens present, 0: brlens absent
#===============================================================================
sub test_has_brlens {

    my ($file) = @_;
    my $brl    = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        chomp;
        if(/^\s*tree/i) {
            if (/\[/) {
                s/(\[.+?\])//g;
            }
            if ( $_ =~ /:/ ) {
                $brl = 1;
            }
            last;
        }
        elsif (/\)\s*;\s*$/) {
            if ( $_ =~ /:/ ) {
                $brl = 1;
            }
            last;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    return $brl;

} # end of test_has_brlens


#===  FUNCTION  ================================================================
#         NAME: test_has_clockrate
#      VERSION: 04/09/2013 01:37:48 PM
#  DESCRIPTION: Tests for presence of clockrate in the tree description.
#               Warning: no error checking. Assumes MrBayes clock tree format.
#               tree gen.5000000[&B Igrbranchlens{all}] = [&R] [&clockrate=8.332129945298162e-04] (1:4.907...
#   PARAMETERS: String containing tree description.
#      RETURNS: 1: clockrate present, 0: clockrate absent
#===============================================================================
sub test_has_clockrate {

    my ($file) = @_;
    my $clr    = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        chomp;
        if(/^\s*tree/i) {
            if ( $_ =~ /&clockrate/ ) {
                $clr = 1;
            }
            last;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    return $clr;

} # end of test_has_clockrate


#===  FUNCTION  ================================================================
#         NAME: test_has_support
#      VERSION: 03/03/2017 12:10:35 PM
#  DESCRIPTION: Tests for presence of support values in the tree description.
#               Warning: no error checking.
#               Assumes Nexus or Newick tree format (tree on single line).
#   PARAMETERS: string containing tree description
#      RETURNS: 1: support values present, 0: support values absent
#===============================================================================
sub test_has_support {

    my ($file)  = @_;
    my $support = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        chomp;
        if((/^\s*tree/i) || (/\)\s*;\s*$/)) {
            if (/\[/) {
                if (/&prob=\d+/i) {
                    $support = 1;
                    last;
                }
                else {
                    s/(\[.+?\])//g;
                }
            }
            if ( $_ =~ /\)\d+/ ) {
                $support = 1;
            }
            last;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    return $support;

} # end of test_has_brlens


#===  FUNCTION  ================================================================
#         NAME: test_figtree_format
#      VERSION: 02/02/2007 12:12:05 AM CET
#  DESCRIPTION: Tests for figtree format, i.e., presence of comments in square
#               brackets within the tree string.
#               Warning: no extensive error checking. Looks only for 'tree' and '&'
#   PARAMETERS: string containing tree description
#      RETURNS: 1: figtree format, 0: not figtree format
#===============================================================================
sub test_figtree_format {

    my ($file)  = @_;
    my $figtree = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        chomp;
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
#      VERSION: Mon 27 maj 2019 14:02:44
#  DESCRIPTION: Documentation
#         TODO: Add examples using rmsupport and about converting .con.tre files
#===============================================================================
POD:
=pod

=head1 NAME

burntrees.pl


=head1 VERSION

Documentation for burntrees.pl version 0.3.1


=head1 SYNOPSIS

burntrees.pl [--burnin=<number>] [--pburnin=<number>] [--start=<number>] [--end=<number>] [--jump=<number>] [--IFeelLucky=<number>] [--treesonly] [--rmbrlens] [--rmcomments] [--rmsupport] [--sci2norm=<nr>] [--seed=<nr>] [--myr] [--[no]close] [--getinfo] [--[no]labels] [--format=altnexus|phylip] [--outfile=<file_name>] FILE [> OUTPUT]


=head1 DESCRIPTION

Script for manipulating tree (*.t, *.trprobs, *.con, *.trees) and parameter (*.p, *.log) files
from MrBayes (v.3), BEAST, and PhyloBayes.

The script extracts trees and (by default) the taxon translation table and the trailing "end;"
from tree file.

A number of options are available:

Any contiguous interval of trees can be printed, as well as trees only (nothing other than tree descriptions).

The samples can be thinned by setting a value for how many trees to jump before next is printed.

Branch lengths, support values, and comments can be removed from trees before printing.

Branch lengths in scientific numeric format can be transformed to a fixed numeric format.

MrBayes clock trees with branch lengths in substitutions per site can be transformed to branch lengths in time units.

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

Print information about the number of trees (or samples in p file), thinning, number of samples, presence of branch lengths and support values in file and quit.


=item B<-h, --help>

Prints help message and exits.


=item B<-i, --IFeelLucky=>I<number>

Specify a probability (value between 0 -- 1) for each tree to be printed. That is, print each tree with prob. I<number>.
The B<-i> option can be combined with the option -B<seed> to create reproducible results.
Note that B<--IFeelLucky> has precedence over B<--jump>.


=item B<-j, --jump=>I<number>

Specify a thinning. That is, print every I<number> tree.


=item B<-l, --labels>

Print trees using sequence (taxon) labels instead of the sequence numbers from the translation table.
B<--nolabels> (which is the default) prevents the sequence numbers to be substituted.


=item B<-ma, --man>

Displays the manual page.


=item B<-my, --myr>

Transform branch lengths in a MrBayes clock tree from substitutions per site to time units.


=item B<-o, --outfile=>I<file_name>

Print directly to file I<file_name> instead of standard out.


=item B<-p, --pburnin=>I<number>

Start printing after a fraction of the run, where I<number> is a percentage (e.g. "50" for half the run).


=item B<-rmb, --rmbrlens>

Remove branch lengths from trees.


=item B<-rmc, --rmcomments>

Remove comments (text within, and the enclosing square brackets) from trees.


=item B<-rms, --rmsupport>

Remove support values (bootstrap/posterior probabilities) from trees.


=item B<-sc, --sci2norm=I<number>>

Translate branch lengths from scientific to normal or fixed. Change the precision by specifying the (optimal) I<number>.


=item B<-se, --seed=I<number>>

Set a seed for the I<-i> option to create reproducible sampling results.


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
  burntrees.pl --rmcomments data.t
  burntrees.pl --ifeellucky=0.50 data.t
  burntrees.pl -b=10 -j=10 -t -rmb data.t
  burntrees.pl --treesonly -b=10 data.p
  burntrees.pl --format=altnexus data.t
  burntrees.pl -f=phylip --outfile=data.phy data.t
  burntrees.pl -f=p -b=1 data.con.tre


For concatenation of several files, use (note the append redirection, ">>". See also separate script catmb.pl):

  burntrees.pl -b=10 -noc  data.run1.t >  data.t
  burntrees.pl -b=40 -t    data.run2.t >> data.t
  burntrees.pl -b=20 -t -c data.run3.t >> data.t

To print the '#NEXUS', 'begin trees;', and the translation table only (no trees), or to print the header in .p, .trace, or .log files, use

  burntrees.pl -i -noc data.t

To convert a MrBayes file to a 'altnexus' file or a 'phylip' file, use

  burntrees.pl --format=altnexus data.t
  burntrees.pl --format=phylip data.t

To concatenate several files in to one altnexus formatted file, use (note the combination of '--format' and '--labels')

  burntrees.pl -b=10 -noc  -f=a data.run1.t >  data.t
  burntrees.pl -b=40 -t    -l   data.run2.t >> data.t
  burntrees.pl -b=20 -t -c -l   data.run3.t >> data.t

To extract the 10th tree in phylip format, use

  burntrees.pl --start=10 --end=10 --format=phylip data.t

To extract the second tree in the MrBayes .con.tre file (simple format) in phylip format, use

  burntrees.pl -b=1 -f=p data.con.tre

To remove comments from tree descriptions (e.g. trees in "figtree" format), use

  burntrees.pl --rmcomments data.con.tre

To remove support values from tree description for the first tree in a MrBayes .con.tre file (simple format) while changng branch lengths from scientific to normal (three decimals), use

  burntrees.pl --start=1 --end=1 --rmsupport --sci2norm=3 data.con.tre

To remove branch lengths from tree descriptions, use

  burntrees.pl --rmbrlens data.t

To change the branch length format from scientific to numerical (three decimals) use

  burntrees.pl --sci2norm=3 data.con.tre

To change the MrBayes clock branch length format from substitutions per site to time units (Myr) use

  burntrees.pl --myr clock.t
  burntrees.pl --myr --rmcomments clock.t

To (randomly) sample a specific number of trees, use burntrees.pl together with 'shuf'. Note that you could sample more trees than you have in your tree file (i.e., performing a bootstrap)!

  burntrees.pl -t data.t | shuf -rn 10
  burntrees.pl -t -f p data.t | shuf -rn 1000

=head1 AUTHOR

Written by Johan A. A. Nylander


=head1 REPORTING BUGS

Please report any bugs to I<Johan.Nylander @ nbis.se>.


=head1 DEPENDENCIES

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


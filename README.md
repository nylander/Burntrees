# Burntrees

Perl script for printing and manipulating output from phylogenetic MCMC
programs.

Currently, the script will handle tree- and parameters files from
[MrBayes](https://github.com/NBISweden/MrBayes) (.t, .con, .p),
[BEAST](http://beast.bio.ed.ac.uk/)/[BEAST2](http://beast2.org/) (.trees,
.log), and [PhyloBayes](http://www.phylobayes.org/) (.trace, .treelist).

The script may also be able to handle other Newick- or Nexus-formatted tree
files.

For usage and examples, please see the Documentation below.


#### Files:

- `burntrees.pl`   -- main script
- `catmb.pl`       -- helper script for concatenating files using
  `burntrees.pl`. In order for `catmb.pl` to work, `burntrees.pl` must be
  installed in the PATH.  Alternatively, the full path to `burntrees.pl` can be
  specified by editing `catmb.pl`.
- `data.t`, `data.p` -- test data


## Documentation


### NAME

burntrees.pl

### VERSION

Documentation for burntrees.pl version 0.3.1

### SYNOPSIS

    burntrees.pl [--burnin=<number>] [--pburnin=<number>] [--start=<number>]
    [--end=<number>] [--jump=<number>] [--IFeelLucky=<number>] [--treesonly]
    [--rmbrlens] [--rmcomments] [--rmsupport] [--sci2norm=<nr>]
    [--seed=<nr>] [--myr] [--[no]close] [--getinfo] [--[no]labels]
    [--format=altnexus|phylip] [--outfile=<file_name>] FILE [> OUTPUT]

### DESCRIPTION

Script for manipulating tree (*.t, *.trprobs, *.con, *.trees) and
parameter (*.p, *.log) files from MrBayes (v.3), BEAST, and PhyloBayes.

The script extracts trees and (by default) the taxon translation table
and the trailing "end;" from tree file.

A number of options are available:

Any contiguous interval of trees can be printed, as well as trees only
(nothing other than tree descriptions).

The samples can be thinned by setting a value for how many trees to jump
before next is printed.

Branch lengths, support values, and comments can be removed from trees
before printing.

Branch lengths in scientific numeric format can be transformed to a
fixed numeric format.

MrBayes clock trees with branch lengths in substitutions per site can be
transformed to branch lengths in time units.

A random set of trees can be printed from the tree file.

Lines can also be extracted from a MrBayes .p file.

Trees can be printed in Phylip (Newick) format or as altnexus (sequence
labels instead of numbers), that is, the script can serve as an
efficient tree format converter.

### OPTIONS

Mandatory arguments to long options are mandatory for short options too

    -b, --burnin=*number*
            Start printing after tree *number*.

    -c, --close
            Forces a trailing "end;" to be printed after the last tree.
            --noclose prevents the "end;" to be printed. Note that the
            trailing "end;" in the tree file is printed by default unless
            --noclose is given.

    -co, --concatenate
            Concatenate several files -- *Not yet implemented*. See USAGE
            for alternatives.

    -e, --end=*number*
            End the printing of trees after tree *number* (inclusively). If
            no --end is given, prints to last tree in file.

    -f, --format=*format*
            Trees are printed as specified by *format*, where *format* is
            either *altnexus*: with sequence (taxon) labels instead of
            numbers, or *phylip* (the Newick format).

    -g, --getinfo
            Print information about the number of trees (or samples in p
            file), thinning, number of samples, presence of branch lengths
            and support values in file and quit.

    -h, --help
            Prints help message and exits.

    -i, --IFeelLucky=*number*
            Specify a probability (value between 0 -- 1) for each tree to be
            printed. That is, print each tree with prob. *number*. The -i
            option can be combined with the option -seed to create
            reproducible results. Note that --IFeelLucky has precedence over
            --jump.

    -j, --jump=*number*
            Specify a thinning. That is, print every *number* tree.

    -l, --labels
            Print trees using sequence (taxon) labels instead of the
            sequence numbers from the translation table. --nolabels (which
            is the default) prevents the sequence numbers to be substituted.

    -ma, --man
            Displays the manual page.

    -my, --myr
            Transform branch lengths in a MrBayes clock tree from
            substitutions per site to time units.

    -o, --outfile=*file_name*
            Print directly to file *file_name* instead of standard out.

    -p, --pburnin=*number*
            Start printing after a fraction of the run, where *number* is a
            percentage (e.g. "50" for half the run).

    -rmb, --rmbrlens
            Remove branch lengths from trees.

    -rmc, --rmcomments
            Remove comments (text within, and the enclosing square brackets)
            from trees.

    -rms, --rmsupport
            Remove support values (bootstrap/posterior probabilities) from
            trees.

    -sc, --sci2norm=*number*
            Translate branch lengths from scientific to normal or fixed.
            Change the precision by specifying the (optimal) *number*.

    -se, --seed=*number*
            Set a seed for the *-i* option to create reproducible sampling
            results.

    -st, --start=*number*
            Start printing from tree *number* (inclusively).

    -t, --treesonly
            Print trees only (do not print taxon descriptions etc.). If used
            on a *.p file it skips the "ID" line and the headers.

    -v, --version
            Prints version message and exits.

    FILE    Reads a Nexus formatted tree FILE, preferrably MrBayes (v.3)
            *.t, *.con, and *.trprobs files. FILE can also be a MrBayes
            parameter file (*.p), or any Nexus formatted tree file, e.g.,
            output from BEAST or other phylogenetic MCMC software (not
            thoroughly tested!).

    OUTPUT  Prints to STDOUT unless --outfile= is used.

### USAGE

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

For concatenation of several files, use (note the append redirection,
">>". See also separate script catmb.pl):

      burntrees.pl -b=10 -noc  data.run1.t >  data.t
      burntrees.pl -b=40 -t    data.run2.t >> data.t
      burntrees.pl -b=20 -t -c data.run3.t >> data.t

To print the '#NEXUS', 'begin trees;', and the translation table only
(no trees), or to print the header in .p, .trace, or .log files, use

      burntrees.pl -i -noc data.t

To convert a MrBayes file to a 'altnexus' file or a 'phylip' file, use

      burntrees.pl --format=altnexus data.t
      burntrees.pl --format=phylip data.t

To concatenate several files in to one altnexus formatted file, use
(note the combination of '--format' and '--labels')

      burntrees.pl -b=10 -noc  -f=a data.run1.t >  data.t
      burntrees.pl -b=40 -t    -l   data.run2.t >> data.t
      burntrees.pl -b=20 -t -c -l   data.run3.t >> data.t

To extract the 10th tree in phylip format, use

      burntrees.pl --start=10 --end=10 --format=phylip data.t

To extract the second tree in the MrBayes .con.tre file (simple format)
in phylip format, use

      burntrees.pl -b=1 -f=p data.con.tre

To remove comments from tree descriptions (e.g. trees in "figtree"
format), use

      burntrees.pl --rmcomments data.con.tre

To remove support values from tree description for the first tree in a
MrBayes .con.tre file (simple format) while changng branch lengths from
scientific to normal (three decimals), use

      burntrees.pl --start=1 --end=1 --rmsupport --sci2norm=3 data.con.tre

To remove branch lengths from tree descriptions, use

      burntrees.pl --rmbrlens data.t

To change the branch length format from scientific to numerical (three
decimals) use

      burntrees.pl --sci2norm=3 data.con.tre

To change the MrBayes clock branch length format from substitutions per
site to time units (Myr) use

      burntrees.pl --myr clock.t
      burntrees.pl --myr --rmcomments clock.t

To (randomly) sample a specific number of trees, use burntrees.pl
together with 'shuf'. Note that you could sample more trees than you
have in your tree file (i.e., performing a bootstrap)!

      burntrees.pl -t data.t | shuf -rn 10
      burntrees.pl -t -f p data.t | shuf -rn 1000

### AUTHOR

Written by Johan A. A. Nylander

### REPORTING BUGS

Please report any bugs to *Johan.Nylander @ nbis.se*.

### DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage

### DOWNLOAD

<https://github.com/nylander/Burntrees>

### LICENSE AND COPYRIGHT

Copyright (c) 2006-2020 Johan Nylander

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

### NAME

catmb.pl

### VERSION

Documentation for catmb.pl version 0.1.2

### SYNOPSIS

    catmb.pl [--burnin=<number>] [--pburnin=<number>] [--jump=<number>]
    FILE FILE [...] [> OUTPUT]

### DESCRIPTION

Script for manipulating tree (*.t) and parameter (*.p) files from
MrBayes (v.3). This script is a helper for concatenating several files
in to one (using the same burnin) using the script 'burntrees.pl'.
burntrees.pl needs to be installed for catmb.pl to work.

### OPTIONS

Mandatory arguments to long options are mandatory for short options too

    -b, --burnin=*number*
            Start printing after tree *number*.

    -p, --pburnin=*number*
            Start printing after a fraction of the run, where *number* is a
            percentage (e.g. "50" for half the run).

    -j, --jump=*number*
            Specify a thinning. That is, print every *number* tree.

    FILE    Reads Nexus formatted tree FILEs, preferrably a MrBayes (v.3)
            *.t file. Can also read and print a MrBayes parameter file (*.p
            file).

    OUTPUT  Prints to STDOUT.

    -h, --help
            Prints help message and exits

    -v, --version
            Prints version message and exits

    -m, --man
            Displays the manual page

### USAGE

Examples:

      catmb.pl --burnin=10 run1.t run2.t run3.t > out.t
      catmb.pl --pburnin=50 run1.t run2.t run3.t > out.t
      catmb.pl --jump=10 run1.t run2.t run3.t > out.t

### AUTHOR

Written by Johan A. A. Nylander

### REPORTING BUGS

Please report any bugs to *Johan.Nylander @ nrm.se*.

### DEPENDENCIES

Needs burntrees.pl to run. Uses Perl modules Getopt::Long and Pod::Usage

### DOWNLOAD

<https://github.com/nylander/Burntrees>

### LICENSE AND COPYRIGHT

Copyright (c) 2006-2020 Johan Nylander

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


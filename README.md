# Burntrees

Perl script for printing and manipulating [MrBayes](https://github.com/NBISweden/MrBayes)
tree files (.t, .con), MrBayes parameter files (.p), or other Nexus-formatted tree files.
Also tested with output from [BEAST](http://beast.bio.ed.ac.uk/)/[BEAST2](http://beast2.org/) (.trees, .log),
and [PhyloBayes](http://www.phylobayes.org/) (.trace, .treelist).

### Written by:

Johan Nylander

### Examples

(See documentaion for more)

- Burnin 10 trees

        burntrees.pl --burnin=10 data.t > out.t

- Burnin 50%

        burntrees.pl --pburnin=50 data.t

- Get trees 11 to 30

        burntrees.pl --start=11 --end=30 data.t

- Thin by sampling every 10th tree

        burntrees.pl --jump=10 data.t

- Get trees only

        burntrees.pl --treesonly data.t

- Get info on samples (after the burnin of 10)

        burntrees.pl --getinfo -b=10 data.t

- Remove branch lengths from trees

        burntrees.pl --rmbrlens data.t

- Remove comments from trees

        burntrees.pl --rmcomments data.t

- Sample every tree with prob 0.5

        burntrees.pl --ifeellucky=0.50 data.t

- Burnin 10, thin by 10, trees only, remove branch lengths

        burntrees.pl -b=10 -j=10 -t -rmb data.t

- Burnin 10 **from a parameter file**

        burntrees.pl -b=10 data.p

- Burnin 10, remove headers **from a parameter file**

        burntrees.pl --treesonly -b=10 data.p

- Convert trees to altnexus format

        burntrees.pl --format=altnexus data.t

- Convert trees to Phylip (Newick) format, save to file

        burntrees.pl -f=phylip --outfile=data.phy data.t

- Get second tree from MrBayes .con.tre file in Phylip format

        burntrees.pl -f=p -b=1 data.con.tre

- Change branch length format from scientific to numerical (three decimals)

        burntrees.pl --sci2norm=3 data.con.tre

### Files:

* `burntrees.pl`   -- the Perl script

* `burntrees.txt`  -- documentation in text format

* `burntrees.html` -- documentation in html format

* `catmb.pl`       -- Perl script for concatenating files using `burntrees.pl`

* `catmb.txt`      -- documentation in text format

* `catmb.html`     -- documentation in html format

* `data.t`, `data.p` -- test data


### Get started

Make sure the file `burntrees.pl` is executable (`chmod +x burntrees.pl`) and

    ./burntrees.pl

---

# CatMB

Perl script for concatenating several MCMC runs.

In order for `catmb.pl` to work, `burntrees.pl` must be installed in the PATH.
Alternatively, the full path to `burntrees.pl` can be specified by editing `catmb.pl`.
(Then `chmod +x catmb.pl`) and

    ./catmb.pl


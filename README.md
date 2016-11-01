# Burntrees

## burntrees.pl:

Perl script for printing and manipulating [MrBayes](https://github.com/NBISweden/MrBayes)
tree files (.t, .con), MrBayes parameter files (.p), or other Nexus-formatted tree files.
Also tested with output from [BEAST](http://beast.bio.ed.ac.uk/)/[BEAST2](http://beast2.org/) (.trees, .log),
and [PhyloBayes](http://www.phylobayes.org/) (.trace, .treelist).

### Written by:

Johan Nylander

### Examples



### Files:

* `burntrees.pl`   -- the Perl script

* `burntrees.txt`  -- documentation in text format

* `burntrees.html` -- documentation in html format

* `catmb.pl`       -- Perl script for concatenating files using `burntrees.pl`

* `catmb.txt`      -- documentation in text format

* `catmb.html`     -- documentation in html format

* `data.t`, `data.p` -- test data

Get started with `burntrees.pl`:

Make sure the file `burntrees.pl` is executable (`chmod +x burntrees.pl`) and

    ./burntrees.pl

---

# CatMB

## catmb.pl:

Perl script for concatenating several MCMC runs.

In order for `catmb.pl` to work, `burntrees.pl` must be installed in the PATH.
Alternatively, the full path to `burntrees.pl` can be specified by editing `catmb.pl`.
(Then `chmod +x catmb.pl`) and

    ./catmb.pl


Burntrees
=========

burntrees.pl:
    Perl script for printing and manipulating MrBayes tree files (*.t, and
    *.con files), MrBayes parameter files (*.p files), or other Nexus-
    formatted tree files. Can also handle output from BEAST (*.trees and *.log files)

Written by:
    Johan Nylander
 
Files:
    burntrees.pl   -- the Perl script

    burntrees.txt  -- documentation in text format

    burntrees.html -- documentation in html format

    catmb.pl       -- Perl script for concatenating files using burntrees.pl

    catmb.txt      -- documentation in text format

    catmb.html     -- documentation in html format

    dat.t, dat.p   -- test data

Get started with burntrees.pl:
    Make sure the file burntrees.pl is executable (chmod +x burntrees.pl) and
    ./burntrees.pl


CatMB
=====

catmb.pl:
    In order for catmb.pl to work, burntrees.pl must be installed in the $PATH.
    Alternatively, the full path to burntrees.pl can be specified by editing catmb.pl.
    (Then chmod +x catmb.pl) and
    ./catmb.pl


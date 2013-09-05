Instructions for setting up Spheral++
=====================================

It is recommended that you set up and build Spheral in a separate directory.
For example:
    cd ..
    mkdir SPHERAL

Spheral is hosted on SourceForge (http://sourceforge.net/projects/spheral/)
and is managed with Mercurial (http://mercurial.selenic.com/). 

Get spheral (read only, recommended):
    hg clone http://hg.code.sf.net/p/spheral/code spheral
Or (read-write):
    hg clone ssh://username@hg.code.sf.net/p/spheral/code spheral

Config spheral:
    cd spheral/src
    ./boot
    ./configure --with-opt=3
Or, with ANEOS support (BYOA):
    cd spheral/src
    ./boot
    ./configure --with-opt=3 --with-aneos --with-aneos-link="-L/path/to/your/ANEOS -lmaneos"

Build spheral:
    cd spheral/src
    make

For more config and build options, see the not-quite-up-to-date config manual:
    cd spheral/doc
    pdflatex Building.tex

Test spheral:
    cd spheral/tests
    ../../bin/ats -e ../../bin/python integration.ats

Update spheral: 
    cd spheral
    hg pull
    hg update

It is recommended that you run your simulations in a separate directory, and keep
the spheral directory clean.

IMPORTANT: Run your scripts using the python interpreter provided by spheral, not
your system's default python. You can set up an alias, e.g.:
    alias spy '/path/to/spheral/bin/python'

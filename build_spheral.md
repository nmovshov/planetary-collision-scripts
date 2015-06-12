Instructions for setting up Spheral++
=====================================

It is recommended that you set up and build Spheral in a separate directory.
For example:

    cd ..
    mkdir SPHERAL
    cd SPHERAL

Spheral is hosted on SourceForge (http://sourceforge.net/projects/spheral/)
and is managed with Mercurial (http://mercurial.selenic.com/). 

Get spheral (read only, recommended):

    hg clone http://hg.code.sf.net/p/spheral/code spheral

Or (read-write):

    hg clone ssh://username@hg.code.sf.net/p/spheral/code spheral

Config spheral:

    cd code/src
    ./boot
    ./configure --prefix=/path/to/SPHERAL --with-opt=2 --without-opensubdiv --with-scipy --with-matplotlib

Or, with ANEOS support (BYOA):

    cd code/src
    ./boot
    ./configure --prefix=/path/to/SPHERAL --with-opt=2 --with-aneos --with-aneos-link="-L/path/to/your/ANEOS -lmaneos"

Build spheral:

    cd code/src
    make

If you are running GNU make (check with `make --version`) on a multi core machine,
you can speed things up a bit:

    cd code/src
    make -j 2
    
For more config and build options, see the not-quite-up-to-date config manual:

    cd code/doc
    pdflatex Building.tex

Test spheral:

    cd code/tests
    ../../bin/ats -e ../../bin/python integration.ats

Update spheral: 

    cd code
    hg pull
    hg update

It is recommended that you run your simulations in a separate directory, and keep
the spheral directory clean.

IMPORTANT: Run your scripts using the python interpreter provided by spheral, not
your system's default python. You can set up an alias, e.g.:

    alias spy '/path/to/SPHERAL/bin/python'

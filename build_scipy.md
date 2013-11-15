Setting up a SciPy stack to work with Spheral's python
------------------------------------------------------

SciPy (http://scipy.org) is a collection of useful enhancements to Python that
makes some common tasks in scientific programming much easier. Some of our
scripts depend on one or more SciPy packages, so if a script fails to run with
an import exception you will need to setup some of the SciPy stack. Note that
even though your system probably already has a version of SciPy installed, you
need to build a separate version that works with the Spheral Python
interpreter instead of the system Python. Here are the packages and how to
install them.

+ NumPy (http://www.numpy.org)
  The fundamental package providing N-dimensional arrays, array manipulation,
  linear algebra, and element-wise math functions, among other things. It is
  included in the Spheral build so you don't need to do anything.

+ SciPy (http://scipy.org/scipylib/)
  Yes it has the same name as the entire stack. This core library provides
  user friendly numerical routines. Install by:

        git clone https://github.com/scipy/scipy.git
        cd scipy
        /path/to/spheral/bin/python setup.py install
  
  and test by starting Spheral's Python and typing `import scipy`.

+ Matplotlib (http://matplotlib.org)
  A 2D plotting library with many routines similar in syntax and usage to
  MATLAB's most common plotting commands. Not as good as the real MATLAB, but
  works in a pinch. Sadly, the github hosted version doesn't work, and
Sourceforge doesn't allow curl downloads. So, you'll need to manually download
a tar.gz file from [the website](http://matplotlib.org/downloads.html). Then:

        tar -xvzf matplotlib-x.x.x.tar.gz
        cd matplotlib-x.x.x
        /path/to/spheral/bin/python setup.py install
        
  and test by starting Spheral's Python and typing `import matplotlib`.

+ IPython (http://ipython.org)
  An enhanced interactive shell. Much better than the default Python
  interactive interpreter. Very useful for debugging scripts. Install by:

        git clone https://github.com/ipython/ipython.git
        cd ipython
        /path/to/spheral/bin/python setup.py install

  This will put an executable `ipython` file in the Spheral bin directory that 
  you can launch from there, or better yer put a convenient shortcut to in your
  rc file:
  
        alias spy "/path/to/spheral/bin/python"
        alias spi "/path/to/spheral/bin/ipython"
        
That's all.

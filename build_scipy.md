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
  included in the default Spheral build, so normally you don't need to do
  anything. But if you see errors relating to BLAS read the troubleshooting
  section below.

+ SciPy (http://scipy.org/scipylib/)
  Yes it has the same name as the entire stack. This core library provides
  user friendly numerical routines. The Spheral build installs it if configured
  with the `--with-scipy` option. Or, you can install without rebuilding Spheral
  by:

        curl -k -L http://sourceforge.net/projects/scipy/files/scipy/0.14.0/scipy-0.14.0.tar.gz/download > scipy-0.14.0.tar.gz
        tar -xvzf scipy-0.14.0.tar.gz
        cd scipy-0.14.0
        /path/to/spheral/bin/python setup.py install
  
  and test by starting Spheral's Python and typing `import scipy`. If you see
  errors relating to BLAS or LAPACK read the troubleshooting section below.

+ Matplotlib (http://matplotlib.org)
  A 2D plotting library with many routines similar in syntax and usage to
  MATLAB's most common plotting commands. Not as good as the real MATLAB, but
  works in a pinch. The Spheral build will install it if configured with the
  `--with-matplotlib` option. To install independently, use:

        curl -k -L http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.3.1/matplotlib-1.3.1.tar.gz/download         > matplotlib-1.3.1.tar.gz
        tar -xvzf matplotlib-1.3.1.tar.gz
        cd matplotlib-1.3.1
        /path/to/spheral/bin/python setup.py install
        
  and test by starting Spheral's Python and typing `import matplotlib`.

+ IPython (http://ipython.org)
  An enhanced interactive shell. Much better than the default Python
  interactive interpreter. Very useful for debugging scripts. Install by:

        git clone https://github.com/ipython/ipython.git
        cd ipython
        git checkout 0.13.x
        /path/to/spheral/bin/python setup.py install

  This will put an executable `ipython` file in the Spheral bin directory that 
  you can launch from there, or better yer put a convenient shortcut to in your
  rc file:
  
        alias spy "/path/to/spheral/bin/python"
        alias spi "/path/to/spheral/bin/ipython"
        
That's all.

### Troubleshooting
+ BLAS and LAPACK
  The NumPy and Scipy libraries count on your system having several static and 
  shared libraries in one of several standard paths. In particular, SciPy requires
  the libraries `libblas.a`, `libblas.so`, `liblapack.a`, `liblapack.so` be available
  while *NumPy* is built.
  
  Typically, these libraries will be found on your system in some standard location, 
  like `/usr/lib` or `/usr/lib64` and the setup scripts will find them easily. However,
  sometimes these libraries are put in a different location, and you need to tell the 
  scripts where to find them. First locate the static libraries by

        locate libblas.a
        locate liblapack.a
        
  and put the locations into the environment variables `BLAS` and `LAPACK`. Then find
  the shared libraries,
  
        locate libblas.so
        locate liblapack.so
        
  and *add* the locations to the environment variable `LD_LIBRARY_PATH`. That should do
  the trick.
  
  
  

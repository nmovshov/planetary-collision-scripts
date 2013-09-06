Instructions for setting up IPython to work with spheral's python
==================================================================

IPython is an enhanced interactive console, or shell, that runs the python
interpreter. It is especially useful when exploring unknown modules, debugging
new modules, or interactive data analysis. There is a very good chance that
your system already has IPython installed, but just like the system installed
python interpreter this IPython shell will not work with the spheral
extensions. Fortunately it is easy to set up an IPython shell that will work
with spheral's modules.

Download IPython:

    git clone https://github.com/ipython/ipython.git

Install using spheral's python:

    cd IPython
    /path/to/spheral/bin/python setup.py install

Set up a convenient shortcut:

    alias spi '/path/to/spheral/bin/ipython'

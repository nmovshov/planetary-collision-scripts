Interactive work and usage demos with Spheral++ objects
=========================================================
It is recommended to run these scripts in an interactive python session.
IPython (http://ipython.org/) is a much better option than the default python
interpreter. You will have to build IPython from source, using Spheral's python
instead of the system python! This is very easy to do.

  - `equations_of_state.py` [BV]
    + To use the ANEOS object remember to build Spheral with ANEOS support
  - `load_interactive_spheral.py`
    + All the components of a spheral run are loaded into the gloabl namespace
    + The Spheral member functions are loaded under the namespace `sph`

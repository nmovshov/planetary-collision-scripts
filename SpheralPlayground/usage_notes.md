Interactive work and usage demos with Spheral++ objects
=========================================================
It is recommended to run these scripts in an interactive python session.
IPython (http://ipython.org/) is a much better option than the default python
interpreter. You will have to build IPython from source, using Spheral's python
instead of the system python! This is very easy to do.

  - `equations_of_state.py` 
    + To use the ANEOS object remember to build Spheral with ANEOS support
    + Currently, ANEOS behaves a little differently from other solid equations
      of state. There is no pressure zeroing mechanism in highly extended
      states, so the way to guard against high negative pressures is to put 
      reasonable limits on density in the interpolation table. Note that this
      will not affect the actual SPH node densities!

  - `load_interactive_spheral.py`
    + All the components of a spheral run are loaded into the gloabl namespace
    + The Spheral member functions are loaded under the namespace `sph`

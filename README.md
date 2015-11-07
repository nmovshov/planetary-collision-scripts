Planetary Collision Scripts
===========================
Naor's Python scripts for initiating, driving, and analyzing spheral simulations.  

**WARNING: work in progress, updated often!**

[Spheral++](http://sourceforge.net/projects/spheral/) is a Lagrangian, ASPH based
hydrocode coupled with an oct-tree gravitational code. Spheral is developed and
maintained by Mike Owen at LLNL. Some reasons that make spheral a good choice for
planetary collision modeling are:  
  - Good scalability (strong and weak) on compute clusters
  - Adaptive sph implementation with ellipsoidal varying smoothing length  
  - Tensor damage model  
  - Multiple equations of state  

This collection of python scripts is an attempt to provide 'templates' to generic
planetary collision simulations including setting up, running, and analyzing
simulation data, so that users can quickly start production runs without the need
to learn the spheral API. For this reason I chose a script-like style of python
modules (where all the action happens in the init section). These scripts are
meant to be used from the shell so they are documented with comments rather than
doc strings. Users typically set a small number of parameters in (a copy of) the
.py files, ideally in no more than the first 20 or so commented lines. I also made
most files executable to simplify the mpirun command. But note that they need to
be executed by the python interpreter built by spheral, not the system python!

To get started read `install.md` and `build_spheral.md`.

Maintained by: Naor Movshovitz (nmovshov@gmail.com)

Directory: **(Work in progress, updated often!)**  
  - `Benchmarks/` (place holder)- Tests of accuracy and performance.  
  - `PlanetFactory/` (place holder) - Scripts for setting up and saving custom
     planets.
  - `GenericCollisions/` - The basic planetary collision scenarios, to be
     customized with problem specific requirements.
  - `Analysis/` - Scripts to load and visualize simulation data.  
  - `SpheralPlayground/` - Various python code chunks for demonstrating usage of
     spheral classes. Usually invoked for interactive work, these should be run
     with iPython if possible.
  - `Help/` - Some usage notes, troubleshooting notes, etc.
  - `shelpers.py` - A module containing some convenience methods used by the other
     scripts. This is a true python module, not an executable.
  - `MATERILAS.md` - Table of material and eos selection for use in many of above
     scripts.
  - `PlanetNodeGenerators.py` - Custom node placement classes. Eventually should
     move into the spheral cose base.

ucsc-spheral-scripts
======================
**WARNING: work in progress, updated often!**

Naor's Python scripts for driving spheral simulations.  

Spheral++ (http://sourceforge.net/projects/spheral/) is a Lagrangian, ASPH based hydrocode coupled with an oct-tree gravitational code. Spheral is developed and maintained by Mike Owen at LLNL. Some reasons that make spheral a good choice for planetary collision modeling are:  
  - Good scalability (strong and weak) on compute clusters
  - Adaptive sph implementation with ellipsoidal varying smoothing length  
  - Tensor damage model  
  - Multiple equations of state  

This collection of python scripts is an attempt to provide 'templates' to generic planetary collision simulations, including setting up, running, and analyzing simulation data, so that users can quickly start production runs without the need to learn the spheral api. For this reason I chose a script-like style of python modules (where all the action happens in the init section). These scripts are meant to be used from the shell so they are documented with comments rather than doc strings. Users typically set a small number of parameters in (a copy of) the .py files, ideally in no more than the first 20 or so commented lines. I also made most files executable to simplify the mpirun command. But note that they need to be executed by the python interpreter built by spheral, not the system python!

Mentained by: Naor Movshovitz (nmovshov@gmail.com)

Directory: **(WARNING: work in progress, updated often!)**  
  [MF=mostly functional; BV=beta version (use with care); ID=in development (do not use); PH=place holder]
  - `Benchmarks` - Tests of accuracy and performance.  
    + `NakamuraFujiwara.py` - Collision of two spheres with strength and damage, no gravity. Try to reproduce Nakamura and Fujiwara (1991). [MF] 
    + `polytropic_planet.py` - Hydrostatic equilibrium of a fluid polytropic planet, compared with known analytic solution. [ID]  
  - `PlanetFactory` - Scripts for setting up and saving sph node lists modeling generic and custom planets.
    + `build_hydrostatic_single_material_planet.py` - Set up a spherical, single material, fluid planet and run to hydrostatic equilibrium. No strength or damage model. [ID]
    + `build_hydrostatic_core_mantle_planet.py` - Set up a two-layer, fluid planet and run to hydrostatic equilibrium. No strength or damage model. [ID]
  - `GenericCollisions` - The basic planetary collision scenarios, to be customized with problem specific requirements.  
    + `strength_collision.py` - Two spheres colliding in the strength regime. Hydro, strength, damage models. No gravity. [BV] 
  - `Analysis` -  Scripts and interactive modules to load and visualize simulation data.  
  - `SpheralPlayground` -Various python code chunks for demonstrating usage of spheral classes. Usually invoked for interactive work, these should be run with iPython if possible.   
    + `gnu_plot.py` - Demonstrating the use the gnuplot python package. [MF]
    + `equations_of_state.py` - Demonstrating the use of spheral's equation-of-state objects. [BV]
    + `load_interactive_spheral.py` - Just the `import` commands that bring spheral classes to the local namespace, for interactive exploration. Works best with iPython. [MF]
  - `shelpers.py` - A module containing some convenience methods used by the other scripts. This is a true python module, not an executable. [ID]
  - `MATERILAS.md` - Table of material and eos selection for use in many of above scripts. [PH]

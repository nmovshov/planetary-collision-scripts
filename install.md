Installation instructions
=========================

There isn't anything to install, but it will help if you add an environment 
variable that holds the absolute path to the base directory (wherever you 
downloaded this package to). This will allow you to run copies of the scripts
from outside the repo clone without modifying anything. Here is my suggested
setup:

+ Use your git client to download the scripts:  
        
        git clone https://github.com/nmovshov/planetary-collision-scripts.git
        
+ We'll call this newly created, top-level directory `<pcs>`. Add the absolute 
  path to this directory to variable called `PCSBASE` in your `.bashrc` or 
  `.cshrc` file:
        
        export PCSBASE = "path/to/<pcs>"

  or  
  
        setenv PCSBASE "path/to/<pcs>"
  
  respectively.
  
+ Now you can run the scripts from where they are located, or copy a script to 
  any location, like `Project_Name/runs`, rename it and edit with your settings,
  and run from there.

+ That's it. But be sure to update the scripts often:

        cd <pcs>
        git pull
  
  and be sure to look in `<pcs>/build_spheral.md` for instructions on
  installing `SPHERAL`.

  In addition, to help with post processing and analysis you may want to integrate
  the SciPy stack with the SPHERAL-built python. Look in
  `<pcs>/Help/build_scipy.md` for instructions.

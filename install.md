Instalation instructions
=========================

There isn't anything to install, but it will help if you add the base directory
to your `PYTHONPATH` environment variable. For completeness, here is what you do
to start using the scripts:

+ Use your git client to download the scripts:  
        
        git clone https://github.com/nmovshov/ucsc-spheral-scripts.git
        
+ We'll call this newly created, top-level directory `<uss>`. Add the path to `<uss>`
  to your `PYTHONPATH` variable, In your `.bashrc` or `.cshrc` file put:
        
        export PYTHONPATH = "path/to/<uss>:$PYTHONPATH"

  or  
  
        setenv PYTHONPATH "path/to/<uss>:$PYTHONPATH"
  
  respectively.
  
+ That's it. But be sure to update the scripts often:

        cd <uss>
        git pull

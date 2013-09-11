Material selection table
==========================
In all scripts provided in this package, the user can choose between the
following equations of state, one per node list, by specifying a material
string.

string          |  EOS used
- - - - - - - - - - - - - - - - - - - - - - - - -
'water'         |  Tillotson parameters for liquid water
'h2oice'        |  Tillotson parameters for solid water ice
'dirtyice'      |  Tillotson parameters for 30% silicate in water ice
'granite'       |  Tillotson parameters for granite
'basalt'        |  Tillotson parameters for basalt
'nylon'         |  Tillotson parameters for nylon (solid)
'SiO2'          |  M/ANEOS SiO2 with Melosh (2007) correction
'poly'          |  A polytrope with n=1 and K=2e5

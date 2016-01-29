Performance scaling of SPHERAL on Hyades
----------------------------------------

Running `<pcs>/GenericCollisions/hydrostatic_single_material_collision.py` on
Hyades with a 1000 km radius ice target and a 500 km radius ice impactor at 3
km/s for 10 cycles.

| # SPH nodes    | 2 cores | 4 cores | 8 cores | 16 cores |
| -------------: | ------: | ------: | ------: | -------: |
|  52480 (nx=40) |   12    |    7    |    5    |    4     |
| 178048 (nx=60) |   

Values are average time interval as reported by the spheral controller at the
end of the run.
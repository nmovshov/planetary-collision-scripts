Performance scaling of SPHERAL on Hyades
----------------------------------------

Running `<pcs>/GenericCollisions/hydrostatic_single_material_collision.py` on
Hyades with a 1000 km radius ice target and a 500 km radius ice impactor at 3
km/s for 10 cycles.

| # SPH nodes (nx, dt)    | 2 cores | 4 cores | 8 cores | 16 cores |
| ---------------------:  | ------: | ------: | ------: | -------: |
|   52480 (nx= 40, dt~4s) |   12    |    7    |    5    |    4     |
|  178048 (nx= 60, dt~3s) |   43    |   25    |   18    |   13     |
|  422980 (nx= 80, dt~2s) |   90    |   49    |   40    |   31     |
|  827630 (nx=100, dt~1s) |  222    |  131    |   87    |   70     |
| 1102095 (nx=110, dt~1s) |  N/A    |  206    |  124    |   88     |
| 1431728 (nx=120, dt~1s) |  418    |  244    |  167    |  114     |

Values are average time interval as reported by the spheral controller at the
end of the run.
Planet building scripts
=========================
Random notes in no particular order.

* If self gravity is dominant you may want to get rid of negative pressures

* When using a stiff eos (so almost always) it's best to update density with
  the continuity equation, i.e. use IntegrateDensity for densityUpdate.

* I suspect an adaptive smoothing length is not really necessary for planet
  building. On the other hand, the plain (non-"A") SPH objects might not have
  all the properties some methods rely on. So the solution is to set the
  minimum allowable hmin/hmax ratio to 1.0.

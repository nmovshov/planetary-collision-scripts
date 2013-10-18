Planet building scripts
=========================
Random notes in no particular order.

* If self gravity is dominant you may want to get rid of negative pressures

* When using a stiff eos (so almost always) it's best to update density with
  the continuity equation, i.e. use IntegrateDensity for densityUpdate.

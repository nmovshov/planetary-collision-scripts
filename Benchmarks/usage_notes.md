Testing stability, accuracy, and performance
==============================================

  - `NakamuraFujiwara.py` 
    + See [Nakamura and Fujiwara, 1991](http://linkinghub.elsevier.com/retrieve/pii/001910359190040Z)

  - `polytropic_planet.py` 
    + For a refresher on polytropes and the solution for a polytropic index n=1
      check out this [summary] (http://www.astro.princeton.edu/~gk/A403/polytrop.pdf).
    + The gist of things is:
      
        ```rho = rho_c * sin(a * r)/(a * r)```
      
      where

        ```a = sqrt(2 * pi * G / K)```
      
      where K is the polytropic constant and G is the gravitational constant. And
      
        ```rho_c = (pi^2 / 3) * rho_av```
         
      where
      
        ```rho_av = 3 * M / (4 * pi * R^3)```
      
      where M is the planet's mass and
      
        ```R = pi / a```
        
      which, you will notice, is independent of mass.

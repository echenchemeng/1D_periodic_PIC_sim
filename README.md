A Julia rewrite of the following project:
https://medium.com/swlh/create-your-own-plasma-pic-simulation-with-python-39145c66578b
https://github.com/pmocz/pic-python


We use a semi-lagrangian method by "shifting" particles around using the classic PIC equations, then rebuild the density using the particles that were shifted around. 
The little velocity perturbation at the beginning gets us a fun little vortex-like pattern quite quickly in the simulation.

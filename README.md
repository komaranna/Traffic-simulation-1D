## Traffic-simulation-1D

/2011-12-12/

C code;

Simulates a 1D traffic situation on a lattice: cars enter on the left with rate alpha, exit to the right with rate beta, and move one to the right with a Poisson rate, i.e. the next time the car moves is taken from an exponential distribution.

N=1000 fixed lattice.

Tlength is the time the simulation stops.

Multi-threaded: nprocs is the number of processors one would like to use.

Outputs 2 text files:

phase_diagram.txt -- average density of cars as a function of alpha and beta;

densityprofileX.txt -- stationary density profile for different alpha and beta values.

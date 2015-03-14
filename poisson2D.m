##	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
##	Students in M1 IHPS 2014 / 2015
##	Project regarding incomplet factorisation (ILU)
##	LANGUAGE USED : Octave

function [A] = poisson2D(nx, ny)
n = nx * ny;
A = 4 * eye(n);
A = A + diag(-1*ones(n-1,1),1) + diag(-1*ones(n-1,1),-1) + diag(-1*ones(n-nx,1),nx) + diag(-1*ones(n-nx,1),-1*nx);

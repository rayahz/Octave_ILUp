##	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
##	Students in M1 IHPS 2014 / 2015
##	Project regarding incomplet factorisation (ILU)
##	LANGUAGE USED : Octave

function [x, iter, r] = CG(A, b)
n = length(b);
x = zeros(n, 1);
iter = 0;
maxiters = max(100, sqrt(n));
r = b';
M = r' * r;
p = r;

while(((norm(r) / norm(b)) > eps) && (iter < maxiters))
	iter = iter + 1;
	alpha = M / (p' * A * p);
	x = x + alpha * p;
	r = r - alpha * A * p;
	Mold = M;
	M = r' * r;
	beta = M / Mold;
	p = r + beta * p;
end;
endfunction

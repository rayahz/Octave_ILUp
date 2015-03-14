##	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
##	Students in M1 IHPS 2014 / 2015
##	Project regarding incomplet factorisation (ILU)
##	LANGUAGE USED : Octave

function [x, iter, r] = PCG(A, b, B)
n = length(b);
x = zeros(n, 1);
iter = 0;
maxiters = max(100, sqrt(n));
r = b';
M = r' * B * r;
p = B * r;

while((norm(r) / norm(b) > eps) && (iter < maxiters))
	iter = iter + 1;
	alpha = M / (p' * A * p);
	x = x + alpha.* p;
	r = r - alpha.* A * p;
	Mold = M;
	M = r' * B * r;
	beta = M / Mold;
	p = B * r + beta.* p;
end;
endfunction

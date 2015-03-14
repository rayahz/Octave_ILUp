##	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
##	Students in M1 IHPS 2014 / 2015
##	Project regarding incomplet factorisation (ILU)
##	LANGUAGE USED : Octave

function [x, iter] = GCR(B, b, x0, max_iter, tolerance)
iter = 0;
j = 1;
erreur = inf;
N = size(B, 1);

if(norm(b) == 0)
	v = zeros(N,1);
	return;
end

r0 = b' - B * x0;
p = zeros(N, max_iter);
Bp = zeros(N, max_iter);
p(:, 1) = r0;

while((erreur >= tolerance) && (iter < max_iter))
	Bp(:, j) = B * p(:, j);
	alpha = ((r0)' * Bp(:, j)) / (Bp(:, j)' * Bp(:, j));
	x0 = x0 + alpha.*p(:, j);
	r0 = r0 - alpha.*Bp(:, j);
	p(:, j + 1) = r0;
	
	for i = 1:1:j
		beta = ((B * r0)' * Bp(:, i)) / (Bp(:, i)' * Bp(:, i));
		p(:, j + 1) = p(:,j + 1) - beta.*p(:, i);
	end;
	
	j = j + 1;
	iter = iter + 1;
	erreur = norm(r0, 2) / norm(b, 2);
end;

x = x0;

endfunction

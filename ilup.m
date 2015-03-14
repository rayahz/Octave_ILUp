##	AUTHORS : Walid HAMMACHE, Yanis KHORSI, Rayhana ZIARA
##	Students in M1 IHPS 2014 / 2015
##	Project regarding incomplet factorisation (ILU)
##	LANGUAGE USED : Octave

function [level,B, L, U] = ilup(A,p)
n = size(A,1);

for i = 1:n, for j = 1:n
	if((abs(A(i,j)) > eps) | (i == j))
		level(i,j) = 0;
	else
		level(i,j) = 1000;
	end
end, end

for i = 2:n
	for k = 1:i-1
		if(level(i,k) <= p)
			A(i,k) = A(i,k) / A(k,k);
			for j = k+1:n
				A(i,j) = A(i,j) - A(i,k) * A(k,j);
				if(abs(A(i,j)) > eps)
					level(i,j) = min(level(i,j), level(i,k) + level(k,j) + 1);
				end
			end
		end
	end

	for j = 1:n
		if(level(i,j) > p)
			A(i,j) = 0;
		end
	end
end

B = A;
L=tril(A,-1+eye(n));
U=triu(A);

endfunction

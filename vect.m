function [b] = vect(n)
for i = 1:n/2
	b(i) = i - 1;
end

if(mod(n,2) == 0)
	for i = n/2:-1:1
		b(n - i+1) = i - 1;
	end
else
	for i = n/2-0.5:-1:0
		b(n - i) = i;
	end
end
endfunction

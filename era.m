function [g,eS] = era (S, it)

% dimension, estimate signal
d = max(size(S));
g = stdnormal_rnd(d,1);
g = g / norm(g);

% iteration
for i=1:it
		G  = fft(g) ;
		eS(i) = sumsq(abs(G)-S);
		G2 = S .* (G ./ abs(G));
		g2 = ifft(G2);
		g  = (g2>=0).*g2;
end

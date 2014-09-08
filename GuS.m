function [g,es,eS] = GuS (s,S,it)

% dimension, estimate signal
d = max(size(s));
g = stdnormal_rnd(d,1) + i * stdnormal_rnd(d,1);
g = g / norm(g);

% iteration
for j=1:it
		G  = fft(g) ;
		eS(j) = sumsq(abs(G)-S);
		G2 = S .* (G ./ abs(G));
		g2 = ifft(G2);
		es(j) = sumsq(abs(g2)-s);
		g  = s .* (g2 ./ abs(g2));
end


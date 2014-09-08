% function [g, es, eS] = GuS (s, S, it)
%
% Using the algorithm of Gerchberg and Saxton to reconstruct a complex signal
% with norm 1.
%
% Input:  's'  is the magnitude of the signal
%         'S'  is the magnitude of the fouriertransformation
%         'it' is the number of iteration before the algorithm terminate
%
% Output: 'g'  is the estimate signal
%         'es' is the quadratic error of the estimate signal in every iteration
%         'eS' is the quadratic error of the fouriertransformation of estimate
%              Signal in every iteration

function [g,es,eS] = GuS (s,S,it)

% dimension, first estimate signal
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


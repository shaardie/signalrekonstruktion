% function [g, eS] = hio (S, it, alpha)
%
% Using the hybrid input output algorithm to reconstruct a real signal with
% norm 1.
%
% Input:  'S'     is the magnitude of the fouriertransformation
%         'it'    is the number of iteration before the algorithm terminate
%         'alpha' is an parameter chosen by the user
%
% Output: 'g'     is the estimate signal
%         'eS'    is the quadratic error of the fouriertransformation of
%                 estimate Signal in every iteration

function [g,eS] = hio (it,S,alpha)

% estimate signal
d = max(size(S));
g = stdnormal_rnd(d,1);
g = g / norm(g);

% iteration
for i=1:it
		G  = fft(g) ;
		eS(i) = sumsq(abs(G)-S);
		G2 = S .* (G ./ abs(G));
		g2 = ifft(G2);
		g  = (g2>=0).* g2 + (g2<0).*(g - alpha^-1 * g2);
end

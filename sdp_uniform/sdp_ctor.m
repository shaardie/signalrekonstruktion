% sdp_ctor construct a real SDP in sedumi-format out of a complex SDP in
% sedumi-format
%
% Input:  complex SDP A,b,c in sedumi-format
% Output: real SDP A,b,c in sedumi-format

function [A,b,c] = sdp_ctor (A,b,c)

   % dimension
   n = sqrt(length(c));
   % number of constraints
   m = size(A,2);
   % new number of constraints
	k = n*(n-1);

   % construct real c out of complex c
   c = reshape(c,n,n);
   c = vec([real(c) , -imag(c) ; imag(c) , real(c)]);

   % construct real b out of complex b
   b = 2 * b;
   b = [b ; zeros(k , 1)];

   % construct real A out of complex A
   % construct real constraints out of complex constraints
   A = [A ; zeros(3*n^2,m)];
   for j = 1:m
      H = reshape(A(1:n^2,j),n,n);
      A(: , j) = vec([real(H) , -imag(H) ; imag(H) , real(H)]);
   end;

   % add control constraints
   A = [A , zeros(4*n^2,k)];
   for i = 1:n-1
      for j = i+1:n
         m++;
         E = (1:n==i).'*(1:n==j) + (1:n==j).'*(1:n==i);
         A(: , m) = vec([E, zeros(n) ; zeros(n) , -E]);
         A(: , m+n*(n-1)/2) = vec([zeros(n), E ; E , zeros(n)]);
      end
   end

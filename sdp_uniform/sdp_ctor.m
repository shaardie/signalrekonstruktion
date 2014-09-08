% sdp_ctor erstellt aus einem komplexen SDP im Sedumi-
% Format ein reelles SDP im Sedumi-Format.
% Input: komplexes SDP A,b,c im Sedumi-Format
% Output: reelles SDP A,b,c im Sedumi-Format
 
function [A,b,c] = sdp_ctor (A,b,c)

   % Dimension 
   n = sqrt(length(c));   
   % Anzahl Constraints
   m = size(A,2);
   % Anzahl der neuen Constraints
   k = n*(n-1);

   % Erstelle aus komplexen c ein reelles c
   c = reshape(c,n,n);
   c = vec([real(c) , -imag(c) ; imag(c) , real(c)]);
   
   % Erstelle aus komplexen b ein reelles b
   b = 2 * b;
   b = [b ; zeros(k , 1)];
   
   % Erstelle aus komplexen A ein reelles A
   % Konstruiere aus den komplexen Constraints
   % reelle Constraints
   A = [A ; zeros(3*n^2,m)];
   for j = 1:m
      H = reshape(A(1:n^2,j),n,n);
      A(: , j) = vec([real(H) , -imag(H) ; imag(H) , real(H)]);
   end;
   
   % Fuege Constraints hinzu, die die Form der Ausgabe
   % kontrolieren
   A = [A , zeros(4*n^2,k)];
   for i = 1:n-1
      for j = i+1:n
         m++;
         E = (1:n==i).'*(1:n==j) + (1:n==j).'*(1:n==i);
         A(: , m) = vec([E, zeros(n) ; zeros(n) , -E]);
         A(: , m+n*(n-1)/2) = vec([zeros(n), E ; E , zeros(n)]);
      end
   end

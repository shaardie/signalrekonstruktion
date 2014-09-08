clear;
% dimension, number of vectors
d = 10;
n = 100;
m = 100;

for l =1:m
	for k =1:n
		% Signal
		s = stdnormal_rnd(d,1) + i * stdnormal_rnd(d,1);
		s = s / norm(s);

		% Erstelle komplexes SDP
		for j = 1:k
   		H(:,j) = stdnormal_rnd(d,1) + i * stdnormal_rnd(d,1) ;
   		H(:,j) = H(:,j) / norm(H(:,j));
			A(:,j) = vec(H(:,j) * H(:,j)');
   		b(j,1) = abs(s' * H(:,j) )^2;
		end

		% Erstelle aus dem komplexen SDP ein reelles SDP
		[A,b,c] = sdp_ctor(A,b,vec(eye(d)));

		% parameter
		K.s = [2*d];
		K.q = [0];
		K.r = [0];
		K.f = [0];

		% Loesen des reellen SDP
		[X,y,z] = csdp(A,b,c,K);

		% Erstelle Matrix mit komplexen Ergebnis
		X = reshape(X,2*d,2*d);
		X = X(1:d,1:d) + i * X(d+1:2*d,1:d);

		% Berechne Fehler
		err(k,l) = norm(X - s*s','fro');

		% Bereinige Daten
		clear H A b K X y z c j s
	end
end

for k=1:n
	meanerr(k) = mean(err(k,:));
end

save sdp_uniform.data

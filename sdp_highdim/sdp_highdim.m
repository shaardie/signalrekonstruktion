clear;
% dimension, number of vectors
d = 10;
n = 60;
m = 100;
max_k = 5;

% parameter
K.s = [d];
K.q = [0];
K.r = [0];
K.f = [0];
c = vec(eye(d));

for k = 2:max_k
	for l =1:m
		for j =10:n
			% Signal
			s = stdnormal_rnd(d,1);
			s = s / norm(s);

			% Erstelle reelles SDP
			for o = 1:j
				H = stdnormal_rnd(d,k);
				H = H / norm(H);
				H = H * (H'*H)^-1 *  H';
				A(:,o) = vec(H);
				b(o,1) = norm(H * s)^2;
			end

			% Loesen des reellen SDP
			[X,y,z] = csdp(A,b,c,K);

			% Erstelle Matrix mit komplexen Ergebnis
			X = reshape(X,d,d);

			% Berechne Fehler
			err(k,j,l) = norm(X - s*s','fro');

			% Bereinige Daten
			clear H A b
		end
	end
end

save sdp_highdim.data

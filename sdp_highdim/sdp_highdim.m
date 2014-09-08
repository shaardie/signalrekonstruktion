% This is a script to reconstruct many real 'd'-dimensional Signals from the
% magnitudes of subspace components for a different numbers of subspaces and
% different dimension of subspaces.

% csdp is needed

clear;
% dimension of the signals
d = 10;
% maximal number of subspaces
n = 60;
% number of reconstructed signals for every number of subspaces
m = 100;
% maximal number of dimensions
max_k = 5;

% parameter in the sedumi-format
K.s = [d];
K.q = [0];
K.r = [0];
K.f = [0];
c = vec(eye(d));

for k = 2:max_k
	for l =1:m
		for j =10:n
			% signal to reconstruct
			s = stdnormal_rnd(d,1);
			s = s / norm(s);

			% construct real SDP
			for o = 1:j
				H = stdnormal_rnd(d,k);
				H = H / norm(H);
				H = H * (H'*H)^-1 *  H';
				A(:,o) = vec(H);
				b(o,1) = norm(H * s)^2;
			end

			% solve real SDP
			[X,y,z] = csdp(A,b,c,K);

			% construct result
			X = reshape(X,d,d);

			% calculate error
			err(k,j,l) = norm(X - s*s','fro');

			% clear data
			clear H A b
		end
	end
end
% save data in file 'spd_highdim.data'
save sdp_highdim.data

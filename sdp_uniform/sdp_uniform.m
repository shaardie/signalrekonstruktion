% This is a script to reconstruct many complex 'd'-dimensional Signals with
% PhaseLift

% csdp is needed

clear;
% dimension of the signals
d = 10;
% maximal number of measurements
n = 100;
% number of reconstructed signals for every number of measurements
m = 100;

for l =1:m
	for k =1:n
		% Signal to reconstruct
		s = stdnormal_rnd(d,1) + i * stdnormal_rnd(d,1);
		s = s / norm(s);

		% construct complex SDP
		for j = 1:k
   		H(:,j) = stdnormal_rnd(d,1) + i * stdnormal_rnd(d,1) ;
   		H(:,j) = H(:,j) / norm(H(:,j));
			A(:,j) = vec(H(:,j) * H(:,j)');
   		b(j,1) = abs(s' * H(:,j) )^2;
		end

		% construct real SDP out of the complex SDP
		[A,b,c] = sdp_ctor(A,b,vec(eye(d)));

		% parameter in sedumi-format
		K.s = [2*d];
		K.q = [0];
		K.r = [0];
		K.f = [0];

		% solve real SDP
		[X,y,z] = csdp(A,b,c,K);

		% construct complex result
		X = reshape(X,2*d,2*d);
		X = X(1:d,1:d) + i * X(d+1:2*d,1:d);

		% calculate error
		err(k,l) = norm(X - s*s','fro');

		% clear data  
		clear H A b K X y z c j s
	end
end

% calulate mean error for every number of measurements
for k=1:n
	meanerr(k) = mean(err(k,:));
end

% save data in file 'spd_uniform.data'
save sdp_uniform.data

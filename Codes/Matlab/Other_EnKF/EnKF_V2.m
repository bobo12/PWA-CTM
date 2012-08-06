function update = EnKF_V2(A, HA, d, R, N, m)

% n = dimension of the state
% m = dimension of the observations
% N = number of ensembles
% A = N samples of the state vector, dim n*N
% H = observation matrix, dim m*n
% d = measurements, dim m*1
% R = covariance matrix for the observations, dim m*m
% HA = H * A; % computes HA, dim m*N
% Ap = A * (eye(N) - ones(N, N) / N); % computes A'

% returns update, dim n*N

HAp = HA * (eye(N) - ones(N, N) / N); % computes HA'
invL = chol(HAp * HAp' + (N - 1) * R, 'lower') \ eye(m); % computes the chol factorization of P
update = A + A * (eye(N) - ones(N, N) / N) * HAp'* (invL' * invL) * ((mvnrnd(d',R,N))' - HA);

end
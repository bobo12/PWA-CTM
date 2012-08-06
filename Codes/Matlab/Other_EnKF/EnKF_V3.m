function update = EnKF_V3(A, HA, d, R, N, m)

% with R diag matrix!

% n = dimension of the state
% m = dimension of the observations
% N = number of ensembles
% A = N samples of the state vector, dim n*N
% H = observation matrix, dim m*n
% d = measurements, dim m*1
% invR = elements of the diagonal of the inverse of the covariance matrix for the observations, dim m*1
% HA = H * A; % computes HA, dim m*N
% Ap = A * (eye(N) - ones(N, N) / N); % computes A'

% returns update, dim n*N
invR = R.^(-1);
HAp = HA * (eye(N) - ones(N, N) / N); % computes HA'
invL = chol((N-1)*eye(N) + HAp' * bsxfun(@times,diag(invR),HAp), 'lower') \ eye(N);
D = (mvnrnd(d',R,N))' - HA;
temp = D - HAp * (invL' * invL) * HAp' * bsxfun(@times,diag(invR),D);
update = A + A * (eye(N) - ones(N, N) / N) * HAp' * bsxfun(@times,diag(invR),temp) / (N-1);
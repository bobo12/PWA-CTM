function  rhoNext = EnKFaprioriMode(rho, percentStateNoise, J, w, d, rhoJ, up, dn)
%
% rho is an ensemble of vectors of the densities on the road at time j
%

% forecast step using godunov scheme
[m,K] = size(rho);% size of the ensemble
% K is the number of ensembles
% taking the mean of rho for the noise
stateNoiseCov = percentStateNoise^2*diag(mean(rho,2).^2); % assuming uncorrelated
stateNoise = mvnrnd(zeros(m,1)',stateNoiseCov,K)'; % state noise ensemble

rhoNext = zeros(m,K);

for l=1:K
    rhoNext(:,l) = godunovSchemeMode(rho(:,l), J, w, d, up, dn)+stateNoise(:,l);
    rhoNext(:,l) = max(min(rhoNext(:,l),rhoJ),0);
end
end
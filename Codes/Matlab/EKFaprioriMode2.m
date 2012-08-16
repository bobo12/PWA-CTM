function  [rhoNext,P] = EKFaprioriMode2(rho, P, percentStateNoise, J, w, d, rhoJ, up, dn)
%
% rho is an ensemble of vectors of the densities on the road at time j
%

% forecast step using godunov scheme
[m,K] = size(rho);% size of the ensemble
% K is the number of ensembles
% taking the mean of rho for the noise
stateNoiseCov = percentStateNoise^2*diag(mean(rho,2).^2); % assuming uncorrelated
stateNoise = mvnrnd(zeros(m,1)',stateNoiseCov,K)'; % state noise ensemble

rhoNext = godunovSchemeMode(rho, J, w, d, up, dn);
rhoNext = max(min(rhoNext,rhoJ),0);

end
function  rhoNext = EnKFaprioriMode(rho,percentStateNoise,dt,dx,v,rhoMax,rhoC)
%
% rho is an ensemble of vectors of the densities on the road at time j
%

% forecast step using godunov scheme
[m,K] = size(rho);% size of the ensemble
% taking the mean of rho for the noise
rhoMean = mean(rho,2);
stateNoiseCov = percentStateNoise^2*diag(rhoMean.^2); % assuming uncorrelated
stateNoise = mvnrnd(zeros(m,1)',stateNoiseCov,K)'; % state noise ensemble

rhoNext = zeros(m,K);
for l=1:K
    rhoNext(:,l) = godunov_scheme(rho(:,l),dt,dx,v,rhoMax,rhoC)+stateNoise(:,l);
    % makes sure no negative value and no value over the jam density
    rhoNext(:,l) = max(rhoNext(:,l),0);
    rhoNext(:,l) = min(rhoNext(:,l),rhoMax);
end
end
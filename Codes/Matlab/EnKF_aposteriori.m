function  rho_next = EnKF_aposteriori(rho,rhoMeasured,H,rhoJam,percentErrorMeasured)

[m,K] = size(rho);
% Compute apriori covariance
covAPriori = cov(rho')';%corresponds to Pj+1
%compute the measurement noise
R = percentErrorMeasured^2*diag(rhoMeasured.^2);
measNoise = mvnrnd(zeros(length(rhoMeasured),1)',R,K)'; % state noise ensemble
% Compute Kalman gain
gain = covAPriori * H' /(H * covAPriori * H' + R);

%Update the forecast : creating aposteriori using measurements
rho_next = rho + gain * (rhoMeasured(:,ones(1,K)) - H * rho + measNoise);%K : ensemble size

rho_next = max(min(rhoJam,rho_next),0);
end
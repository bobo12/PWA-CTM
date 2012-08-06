function ens_rhoi_aposteriori = EnKF_V1(ens_rhoi_apriori, mes_rhoi, cov_mes_rhoi, H)
% K is the number of ensemble realizations
% Computing the ensemble covariance
covr_ens_rhoi_apriori = cov(ens_rhoi_apriori);

% Analysis ::
    % Compute the measurement as an ensemble (of 100)
    ens_mes_rhoi = mvnrnd(mes_rhoi,cov_mes_rhoi,100);
    % Compute Kalman Gain thanks to to measurements 
    G_ens = covr_ens_rhoi_apriori*H'*pinv(H*covr_ens_rhoi_apriori*H' + cov_mes_rhoi);
    % Update the a-priori value to an a-posteriori value 
    ksi_ens = % white zero-mean observation noise with covariance cov_mes_rhoi
    ens_rhoi_aposteriori = ens_rhoi_apriori + G_ens*(ens_mes_rhoi - H*ens_rhoi_apriori + ksi_ens);
end

    
 




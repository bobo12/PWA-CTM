function  rho_next = EKFaposterioriMode(rho,P,rhoMeasured,H,rhoJam,percentErrorMeasured)


%compute the measurement noise
R = percentErrorMeasured^2*diag(rhoMeasured.^2);
% Compute Kalman gain
gain = P * H' /(H * P * H' + R);

%Update the forecast : creating aposteriori using measurements
rho_next = rho + gain * (rhoMeasured - H * rho);
rho_next = max(min(rhoJam,rho_next),0);
end
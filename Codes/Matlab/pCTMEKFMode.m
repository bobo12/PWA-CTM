function rho = pCTMEKFMode(rho_initial, nbTimeSteps, dt, route, ...
    vff, rhoJ, rhoC)

rho = zeros(length(rho_initial),nbTimeSteps);
rhoEKF = rho_initial;
percent_dev_state = 0.01;
P = 0.05^2 * eye(length(rho_initial));

percent_dev_meas = 100; % percentage of deviation of the measurements values
counterEnkfSteps = 0;

alpha = dt / mean(route.cellLength);
omega_f = vff * rhoC / (rhoJ - rhoC);

J = [[0, 1-alpha*omega_f, alpha*omega_f];
    [0, 1-alpha*omega_f, 0];
    [0, 1, alpha*omega_f];
    [0, 1-alpha*vff, 0];
    [alpha*vff, 1, alpha*omega_f];
    [alpha*vff, 1, 0];
    [alpha*vff, 1-alpha*vff, 0]];

w = [0;
    alpha*omega_f*rhoC;
    -alpha*omega_f*rhoC;
    alpha*vff*rhoC;
    -alpha*omega_f*rhoJ;
    -alpha*vff*rhoC;
    0];

d = [[(rhoJ - rhoC) / rhoC, 1, -rhoJ];
    [1, 0, -rhoC];
    [0, 1, -rhoC]];

nbTimeSteps

for j=1:nbTimeSteps
    j
    %Forecast step
    [rhoEKF, P] = EKFaprioriMode(rhoEKF, P, percent_dev_state, J, w, d, rhoJ);% a priori
    if(mod(j-1,6)==0)
        %     if(mod(j-1,12)==0)
        counterEnkfSteps = counterEnkfSteps+1;
        %disp(counterEnkfSteps);
        %get H
        nbRows = size(route.activeSensors{counterEnkfSteps},1);
        if(nbRows~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{counterEnkfSteps},:);
            Hj = [zeros(size(Hj,1),1), Hj, zeros(size(Hj,1),1)];
            rhoEKF = EKFaposterioriMode(rhoEKF,P,Hj*[0;route.densityMeasured(:,counterEnkfSteps);0],Hj,rhoJ,percent_dev_meas);% a posteriori
            %             display('mean_rho');
            %             disp(mean(ensRhoJ,2));
            %store rho
        end
    end
    rho(:,j) = rhoEKF;
end

end
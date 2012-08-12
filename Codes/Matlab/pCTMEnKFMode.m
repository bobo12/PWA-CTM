function rho = pCTMEnKFMode(rho_initial, nbEnsembles, nbTimeSteps, dt, route, ...
            vff, rhoJ, rhoC, useEfficientEnKF, up, dn)

% Generate the initial ensemble
percent_dev_ens = 0.05; % percentage of deviation of the initial values
ens_rho_initial = ...
    mvnrnd(rho_initial',...
           percent_dev_ens^2 * eye(length(rho_initial)),...
           nbEnsembles)';...
%assuming a 5% error on the knowledge of the initial conditon

%Initialize
ensRhoJ = ens_rho_initial;
%Simulate state noise
% Noises
percent_dev_meas = 0.05; % percentage of deviation of the measurements values
percent_dev_state = 0.01; % percentage of deviation of the state values

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

f = [0;
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
    ensRhoJ = EnKFaprioriMode(ensRhoJ, percent_dev_state, J, f, d, rhoJ, up(j), dn(j));% a priori
    if(mod(j-1,6)==0)
%     if(mod(j-1,12)==0)
        counterEnkfSteps = counterEnkfSteps+1;
        %disp(counterEnkfSteps);
        %get H
        nbRows = size(route.activeSensors{counterEnkfSteps},1);
        if(nbRows~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{counterEnkfSteps},:);
            %             update = EnKF_V3(A, HA, invR, N, m)
            if(useEfficientEnKF)
                HA = Hj*ensRhoJ;
                Rj = 0.0009*diag(route.densityMeasured(route.sensorCellMap(...
                    route.activeSensors{counterEnkfSteps}),counterEnkfSteps).^2);
                ensRhoJ = EnKF_V2(ensRhoJ,HA,...
                                    route.densityMeasured(...
                                        route.sensorCellMap(...
                    route.activeSensors{counterEnkfSteps}),counterEnkfSteps),...
                    Rj,nbEnsembles,nbRows);           
            else
                %         %Obtain measurements
                %         rho_mes_j = route.densityMeasured(:,counterEnkfSteps)
                %measurements noise
                % ens_eps = mvnrnd(zeros(1,nbRows),Rj,nbEnsembles)';
                %         without noise ens_eps = 0*mvnrnd(zeros(1,nbRows),Rj,nbEnsembles)';
                %Update using measurements
                ensRhoJ = EnKF_aposteriori(ensRhoJ,Hj*route.densityMeasured(:,counterEnkfSteps),Hj,rhoJ,percent_dev_meas);% a posteriori
                %             display('mean_rho');
                %             disp(mean(ensRhoJ,2));
                %store rho
            end
        end
    end
    rho(:,j) = mean(ensRhoJ,2);
end
function rho = pCTMEnKF(rho_initial, nbEnsembles, nbTimeSteps, dt, route, ...
            vff, rho_max, rho_c, useEfficientEnKF)

% Generate the initial ensemble
percent_dev_ens = 0.05; % percentage of deviation of the initial values
ens_rho_initial = ...
    mvnrnd(rho_initial',...
           percent_dev_ens^2 * eye(length(rho_initial)),...
           nbEnsembles)';...
%assuming a 5% error on the knowledge of the initial conditon

%Initialize
ens_rho_j = ens_rho_initial;
%Simulate state noise
% Noises
percent_dev_meas = 0.05; % percentage of deviation of the measurements values
percent_dev_state = 0.01; % percentage of deviation of the state values

counterEnkfSteps = 0;

for j=1:nbTimeSteps
    %Forecast step
    ens_rho_j = EnKF_apriori(ens_rho_j,percent_dev_state,dt,route.cellLength',vff,rho_max,rho_c);% a priori
    if(mod(j-1,6)==0)
%     if(mod(j-1,12)==0)
        counterEnkfSteps = counterEnkfSteps+1;
        disp(counterEnkfSteps);
        %get H
        nbRows = size(route.activeSensors{counterEnkfSteps},1);
        if(nbRows~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{counterEnkfSteps},:);
            Hj = [zeros(size(Hj,1),1), Hj, zeros(size(Hj,1),1)];
            %             update = EnKF_V3(A, HA, invR, N, m)
            if(useEfficientEnKF)
                HA = Hj*ens_rho_j;
                Rj = 0.0009*diag(route.densityMeasured(route.sensorCellMap(...
                    route.activeSensors{counterEnkfSteps}),counterEnkfSteps).^2);
                ens_rho_j = EnKF_V2(ens_rho_j,HA,...
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
                ens_rho_j = EnKF_aposteriori(ens_rho_j,Hj*[0;route.densityMeasured(:,counterEnkfSteps);0],Hj,rho_max,percent_dev_meas);% a posteriori
                %             display('mean_rho');
                %             disp(mean(ens_rho_j,2));
                %store rho
            end
        end
    end
    rho(:,j) = mean(ens_rho_j,2);
end
function [rho,vel] = pCTM(route, vff, rhoJ, rhoC, dt, ...
                          nbEnsembles, useEfficientEnKF, plotFig, algorithm)

% p-CTM on a fixed route for one day of Pems data

if nargin == 0
    route = buildRoute2;
    %get fundamental diagrams by defining rho_max and v
    vff=28.13;%cf gunes
    rhoJ=1/7;%cf gunes

    %hardcoded (it satisfies the CLF condition for our route    
    useEfficientEnKF = false;
    dt = 5;
    nbEnsembles = 100;
    plotFig = true;
end

tic
nbTimeSteps = floor(route.totalSec/dt);

%initialize rho    
rho_initial = 0.01*ones(route.nbCells + 2, 1);
up = zeros(1,nbTimeSteps);
dn = zeros(1,nbTimeSteps);
%dn(floor(nbTimeSteps/2):nbTimeSteps) = rhoJ * ones(nbTimeSteps-floor(nbTimeSteps/2)+1,1);

%% solving using the godunov scheme :: No Kalman Filter
%

%% solving using godunov scheme and Kalman filter to update

%Assuming the measurements are given by the analytical solution, we
%build the 'measures' vector

rho = zeros(length(rho_initial),nbTimeSteps);
rho(:,1) = rho_initial;
rho(1,:) = up;
rho(length(rho_initial),:) = dn;

switch algorithm
    case 'EnKF'
        rho = pCTMEnKF(rho_initial, nbEnsembles, nbTimeSteps, dt, route, ...
            vff, rhoJ, rhoC, useEfficientEnKF);
    case 'EnKFmode'
        rho = pCTMEnKFMode(rho_initial, nbEnsembles, nbTimeSteps, dt, route, ...
            vff, rhoJ, rhoC, useEfficientEnKF);
    case 'EKFmode'
        rho = pCTMEKFMode(rho_initial, nbTimeSteps, dt, route, ...
            vff, rhoJ, rhoC);
        
end

if plotFig
    figure(1)
    surf(rho(2:end-1,:),'Linestyle','None');
    view(2)
    figure(2)
    vel = vff*(1-rho/rhoJ);
    vel(vel<1e-2)=0;
    surf(vel(2:end-1,:),'Linestyle','None');
    view(2)
    cmap = colormap;
    cmap = flipud(cmap);
    colormap(cmap);
end

toc
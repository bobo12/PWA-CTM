% To use this code, establish an ssh tunnel to mmdevdb1 through
% mmwiki.ccit.berkeley.edu or mmdev.ccit.berkeley.edu

clear all
close all
clc


% Fundamental Diagram parameters
RHO_MAX = 1/7;
RHO_C = RHO_MAX/3;
VFF = 28.13;

% Model timestep
DT = 5;
% EnKF number of ensembles
NUM_ENSEMBLES = 100;
% Use efficient EnKF?
USE_EFFICIENT_ENKF = false;

% Plot density?
PLOT_FIG = true;

ALGORITHM = 'EnKFmode';

route = struct();
route.cellLength = 200 * ones(1,100);
route.nbCells = 100;
route.activeSensors = cell(1000,1);
route.totalSec = 1800;
route.observationMatrix = [];
route.densityMeasured = [];
route.sensorCellMap = [];

[rho,vel] = pCTM(route, VFF, RHO_MAX, RHO_C, DT, NUM_ENSEMBLES, ...
                 USE_EFFICIENT_ENKF, PLOT_FIG, ALGORITHM);
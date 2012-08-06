% To use this code, establish an ssh tunnel to mmdevdb1 through
% mmwiki.ccit.berkeley.edu or mmdev.ccit.berkeley.edu

clear all
close all
clc

% parameters
% Start date
%I880 Bluetooth first deployment
START_TIME = '2012-03-05 00:00:00.000';
% % I15 Bluetooth second deployment
% START_TIME = '2012-04-04 00:00:00.000';
% End date

%I880 Bluetooth first deployment
END_TIME = '2012-03-06 00:00:00.000';
% % I15 Bluetooth second deployment
% END_TIME = '2012-04-05 00:00:00.000';
% Every day will "start" at this hour (anything before this hour
% will be ignored)

START_HOUR = 7;
% Every day will "end" at this hour (anything after this hour will
% be ignored - NOTE: the hour is inclusive, e.g. in this case, the
% last time sample considered is 20:59:30)
END_HOUR = 20;

% Direction of the route
%I880 Bluetooth first deployment
DIRECTION = 'N';
% %I15 Bluetooth first deployment
% DIRECTION = 'N';

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

% MG link ids have to be in the same geographical order as the
% route - for now, the links are chosen manually

% I880 Bluetooth first deployment
ROUTE_LINK_IDS = [186836 315387 1238 150305 233028 97544 47187 ...
                  129727 191223 191222 47190 259805 97540 190831 ...
                  190830 41558 320521 320522 320523 150310 150311 ...
                  129730 83134 151055 101004 278401 150937 150944 155637 ...
                  259138];%ids of model graph links

% %I15 Bluetooth second deployment
% ROUTE_LINK_IDS = [290983,272019,272009,...
%     19940,12032,176773,215609,153957,153949,317341,108776,...
%     129173,153954,215646,215751,150123,225079,319695,319694,...
%     99794,236088,153951,153945,84412,3315,84421];
           
              
% Network ID
%I880 Bluetooth first deployment
NID = 132;
% % I15 Bluetooth second deployment
% NID = 202;

% Plot density?
PLOT_FIG = true;

ALGORITHM = 'EnKF';

route = buildRoute(ROUTE_LINK_IDS, NID, START_TIME, END_TIME, ...
                    DIRECTION, RHO_MAX, START_HOUR, END_HOUR);

[rho,vel] = pCTM(route, VFF, RHO_MAX, RHO_C, DT, NUM_ENSEMBLES, ...
                 USE_EFFICIENT_ENKF, PLOT_FIG, ALGORITHM);
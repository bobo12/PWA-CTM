clc
clear
tic
% This script is aimed at checking the godunov scheme function on a know
% example. CF EE291 lecture

%Global variables
L = 35; 
T =100; 
h = L/350; 
x = (-10:h:L-10)';
dx = h*ones(length(x),1);
rhomax = 4; % 
v = 1; % 
% dx fixed, we defined dt to respect the CFL condition
dt = 0.8*min(dx)/v;
r = dt./dx;

% Boundary & Initial conditions :
% We start from zeros and see the convergence of the EnKF
rho_a = zeros(1,floor(T/dt)); % left phantom-cell
rho_b = zeros(1,floor(T/dt)); % right phantom-cell
rho_initial = 2*ones(length(dx),1);

% initializing rho
rho = zeros(length(rho_initial),floor(T/dt));
% modulable in cas we want to help the EnKF we better initialization
rho(:,1) = rho_initial;
rho(1,:) = rho_a;
rho(length(rho_initial),:) = rho_b;

%% Create the Analytical solution 
% Key points in the time-space diagram to build our analytical solution
xd = (40*v-14)/9;
td = (194-40*v)/(9*v);
xb = L-10;
tb = 2*(xb-20)/v;
% Define rho_analytical depending on the zone of (x,t) we are in
rho_analytical = zeros(length(dx),floor(T/dt));
for i=1:length(dx)
    k = i-10/h+1;
    for j=1:floor(T/dt)
        if(j*dt<=2/v*(k*h-20))
            rho_analytical(i,j) = rhomax/4;
        else if(j*dt>=-1/v*(k*h-20)&& k*h>=20+1/2*j*dt-sqrt(77*j*dt/2))
                rho_analytical(i,j) = 2*(1-(k*h-20)/(v*j*dt));
            else if(j*dt>= 40/(3*v)-4/v*(k*h-23/6) && j*dt>=1/v-2/v*(k*h-10) && j*dt>=-1/v*(k*h-11) && j*dt<td)
                    rho_analytical(i,j) = rhomax;
                else if(k*h>=10)
                        rho_analytical(i,j)= 2*(1+(k*h-10)/(1-v*dt*j));
                    else if(j*dt<=2/v+4/v*(k*h-1) && k*h>=1)
                            rho_analytical(i,j) = rhomax/2;
                        else if(j*dt<=2/v*k*h)
                                rho_analytical(i,j)=(1+k*h-v*j*dt)/(1-v*j*dt/2);
                            else
                                rho_analytical(i,j)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Defines the measurements
% We use the values of the analytical solution as the measurements
% detect_percent : percentage of the cell at which we measure density. We
% assume we don't have loop-detector at the boundaries
% time_percent : percentage of time steps at which we get measurements
% 
%
detect_percent = 0.1; % one detector every 10 cells
time_percent = 0.1; %one measurement every 5 time step


rho_meas = zeros(floor((length(dx)-2)*detect_percent),floor(T/dt*time_percent));
for i=1:floor((length(dx)-2)*detect_percent)
    k = floor(1/detect_percent*i-10/h)+1;
    for j=1:floor(T/dt*time_percent)+1
        l = floor(1/time_percent*j)+1;
        if(l*dt<=2/v*(k*h-20))
            rho_meas(i,j) = rhomax/4;
        else if(l*dt>=-1/v*(k*h-20)&& k*h>=20+1/2*l*dt-sqrt(77*l*dt/2))
                rho_meas(i,j) = 2*(1-(k*h-20)/(v*l*dt));
            else if(l*dt>= 40/(3*v)-4/v*(k*h-23/6) && l*dt>=1/v-2/v*(k*h-10) && l*dt>=-1/v*(k*h-11) && l*dt<td)
                    rho_meas(i,j) = rhomax;
                else if(k*h>=10)
                        rho_meas(i,j)= 2*(1+(k*h-10)/(1-v*dt*l));
                    else if(l*dt<=2/v+4/v*(k*h-1) && k*h>=1)
                            rho_meas(i,j) = rhomax/2;
                        else if(l*dt<=2/v*k*h)
                                rho_meas(i,j)=(1+k*h-v*l*dt)/(1-v*l*dt/2);
                            else
                                rho_meas(i,j)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Building the H matrix : for observations : here it is constant, in the
%% example with real data it will be time-varying

% sparse structure, only ones where there is a detector.
H = zeros(floor((length(dx)-2)*detect_percent),length(dx));
for i=1:floor((length(dx)-2)*detect_percent)
    H(i,floor(1/detect_percent*i)+1)=1;
end

    

%% ENSEMBLE KALMAN FILTER ALGORITHM
K = 100; % ensemble size
% Generate the initial ensemble
percent_dev_ens = 0.01; % percentage of deviation of the initial values
ens_rho_initial = mvnrnd(rho_initial',percent_dev_ens^2*eye(length(rho_initial)),K)';
% Initialize
ens_rho = ens_rho_initial; % ensemble of densities that we will propagate
% Noises
percent_dev_meas = 0.05; % percentage of deviation of the measurements values
percent_dev_state = 0.01; % percentage of deviation of the state values
counter = 1;
for j=1:floor(T/dt)-1  
    % Forecast step : apriori estimate
    ens_rho = EnKF_apriori(ens_rho,percent_dev_state,dt,dx,v,rhomax,'Greenshields');
    % Test if we have a measurement at this time step , in this case
    % ANALYSIS STEP of ENKF
    if(mod(j-1,floor(1/time_percent))==0)
        rho_meas_j = rho_meas(:,counter);
        counter=counter+1;
        % Update using measurements
        ens_rho = EnKF_aposteriori(ens_rho,rho_meas_j,H,rhomax,percent_dev_meas);
    end
        % store rho
    rho(:,j+1) = mean(ens_rho,2);
end
figure('Name','Fixed_steps_and_noise_meas_in_inN_01')
surf(rho'-rho_analytical','Linestyle','None')
view(2)
figure('Name','Fixed_steps_and_noise_meas_in_inN_01_3D')
surf(rho','Linestyle','None')

percent_dev_ens = 0.05;
ens_rho_initial = mvnrnd(rho_initial',percent_dev_ens^2*eye(length(rho_initial)),K)';
% Initialize
ens_rho = ens_rho_initial; % ensemble of densities that we will propagate

counter = 1;
for j=1:floor(T/dt)-1  
    % Forecast step : apriori estimate
    ens_rho = EnKF_apriori(ens_rho,percent_dev_state,dt,dx,v,rhomax,'Greenshields');
    % Test if we have a measurement at this time step , in this case
    % ANALYSIS STEP of ENKF
    if(mod(j-1,floor(1/time_percent))==0)
        rho_meas_j = rho_meas(:,counter);
        counter=counter+1;
        % Update using measurements
        ens_rho = EnKF_aposteriori(ens_rho,rho_meas_j,H,rhomax,percent_dev_meas);
    end
        % store rho
    rho(:,j+1) = mean(ens_rho,2);
end
figure('Name','Fixed_steps_and_noise_meas_in_inN_05')
surf(rho'-rho_analytical','Linestyle','None')
view(2)
figure('Name','Fixed_steps_and_noise_meas_in_inN_05_3D')
surf(rho','Linestyle','None')

percent_dev_ens = 0.1;
ens_rho_initial = mvnrnd(rho_initial',percent_dev_ens^2*eye(length(rho_initial)),K)';
% Initialize
ens_rho = ens_rho_initial; % ensemble of densities that we will propagate

counter = 1;
for j=1:floor(T/dt)-1  
    % Forecast step : apriori estimate
    ens_rho = EnKF_apriori(ens_rho,percent_dev_state,dt,dx,v,rhomax,'Greenshields');
    % Test if we have a measurement at this time step , in this case
    % ANALYSIS STEP of ENKF
    if(mod(j-1,floor(1/time_percent))==0)
        rho_meas_j = rho_meas(:,counter);
        counter=counter+1;
        % Update using measurements
        ens_rho = EnKF_aposteriori(ens_rho,rho_meas_j,H,rhomax,percent_dev_meas);
    end
        % store rho
    rho(:,j+1) = mean(ens_rho,2);
end
figure('Name','Fixed_steps_and_noise_meas_in_inN_1')
surf(rho'-rho_analytical','Linestyle','None')
view(2)
figure('Name','Fixed_steps_and_noise_meas_in_inN_1_3D')
surf(rho','Linestyle','None')

for i = 1:6
        f=figure(i);
        colorbar 
        saveName = get(f,'Name');
        saveas(gcf, ['..\..\Simulation_Results\Test_EnKF\' saveName ], 'png');
        saveas(gcf, ['..\..\Simulation_Results\Test_EnKF\' saveName ], 'fig');
        close
end
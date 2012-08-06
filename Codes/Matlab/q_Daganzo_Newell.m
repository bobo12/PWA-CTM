function qDN_rho = q_Daganzo_Newell(rho,vf,rho_jam,rho_c)
% 
% The inputs of this function are :
%     a vector of densities rho
%     the free-flow velocity vf,
%     the jam density rho_jam
%     the critical density - that separates the flow into two regimes : 
%       free and congestion - rho_c and the density at time t and location
%       x (scalar) rho_x_t
% 
% The outputs is a vector of flows
%
qDN_rho = zeros(length(rho),1);
wf = vf/(rho_jam/rho_c-1); % backwards propagating wave speed
for i=1:length(rho)
    if(rho(i)<0 || rho(i)>rho_jam)
    else if(rho(i)>rho_c)
        qDN_rho(i) = -wf*(rho(i) - rho_jam);
        else
            qDN_rho(i) = vf*rho(i);
        end
    end
end

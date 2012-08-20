function rho_next = godunov_scheme(rho,dt,dx,vf,rhoJam,rhoC)
%depending on the number of inputs in the function
%note that dx is a vector containing the ith cell's length

%r = dt./[dx(1);dx;dx(end)];
r = dt/mean(dx);
%size(rho)
%size(r)
%size(q_Godunov(rho,vf,rhoJam,rhoC))

%rho_next = rho - r.*q_Godunov(rho,vf,rhoJam,rhoC);
rho_next = rho - r * q_Godunov(rho,vf,rhoJam,rhoC);
%TODO: reecrire le code pour rho_next without using q_Godunov

rho_next = max(min(rhoJam,rho_next),0);% just because we don't want negative values
end



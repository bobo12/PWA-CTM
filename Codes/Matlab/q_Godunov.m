function qG = q_Godunov(rho,vf,rhoJ,rhoC)
% function that calculates the modified fundamental diagram while solving
% using the godunov scheme

qG_i_nexti = zeros(length(rho),1);
qG_previ_i = zeros(length(rho),1);

q_rhoC = q_Daganzo_Newell(rhoC,vf,rhoJ,rhoC);
q_rho = q_Daganzo_Newell(rho,vf,rhoJ,rhoC);
       
for i=2:length(rho)-1
    qG_i_nexti(i) = compare_rho(rho(i),rho(i+1),rhoC,q_rho(i),q_rho(i+1),q_rhoC);
    qG_previ_i(i) = compare_rho(rho(i-1),rho(i),rhoC,q_rho(i-1),q_rho(i),q_rhoC);
end            
    qG = qG_i_nexti - qG_previ_i;
end


   
        
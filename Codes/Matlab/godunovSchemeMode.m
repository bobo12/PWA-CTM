function rho_next = godunovSchemeMode(rho, J, f, d, up, dn)
%depending on the number of inputs in the function
%note that dx is a vector containing the ith cell's length

nExt = length(rho);
s = zeros(nExt-1, 1);
m = zeros(nExt-2, 1);
rho_next = zeros(nExt, 1);
rho_next(1) = up;
rho_next(nExt) = dn;

for i=1:nExt-1
    ind = d * [rho(i); rho(i+1); 1] > 0;
    s(i) = find([ind(1)*ind(3); ind(2)*(1-ind(3)); (1-ind(1))*(1-ind(2))]);
end

for i=1:nExt-2    
    m(i) = find([(s(i)==1 && s(i+1)==1);
        (s(i)==1 && s(i+1)==2);
        (s(i)==2 && s(i+1)==1);
        (s(i)==2 && s(i+1)==3);
        (s(i)==3 && s(i+1)==1);
        (s(i)==3 && s(i+1)==2);
        (s(i)==3 && s(i+1)==3)]);
end

for i=2:nExt-1
    rho_next(i) = J(m(i-1),:) * [rho(i-1); rho(i); rho(i+1)] + f(m(i-1));
end

end

% function qG = q_Godunov(rho,vf,rhoJam,rhoC)
% % function that calculates the modified fundamental diagram while solving
% % using the godunov scheme
% 
% qG_i_nexti = zeros(length(rho),1);
% qG_previ_i = zeros(length(rho),1);
% 
% q_rhoC = q_Daganzo_Newell(rhoC,vf,rhoJam,rhoC);
% q_rho = q_Daganzo_Newell(rho,vf,rhoJam,rhoC);
%        
% for i=2:length(rho)+1
%     qG_i_nexti(i) = compare_rho(rho(i),rho(i+1),rhoC,q_rho(i),q_rho(i+1),q_rhoC);
%     qG_previ_i(i) = compare_rho(rho(i-1),rho(i),rhoC,q_rho(i-1),q_rho(i),q_rhoC);
% end            
%     qG = qG_i_nexti - qG_previ_i;
% end
% 
% 
%    
%         
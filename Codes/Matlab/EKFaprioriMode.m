function  [rhoNext, Pnext] = EKFaprioriMode(rho, P, percentStateNoise, J, w, d, rhoJ)
%
% rho is an ensemble of vectors of the densities on the road at time j
%
nExt = length(rho);
% forecast step using godunov scheme
[rhoNext, m] = godunovSchemeMode(rho, J, w, d);

Pnext = zeros(nExt, nExt);

rhoNext = max(min(rhoNext,rhoJ),0);

for i=2:nExt-1
    for j=2:nExt-1
        Pnext(i,j) = J(m(i-1),:) * [P(i-1,j); P(i,j); P(i+1,j)];
    end
end

for i=2:nExt-1
    for j=2:nExt-1
        Pnext(i,j) = [Pnext(i,j-1), Pnext(i,j), Pnext(i,j+1)] * J(m(j-1),:)';
    end
end

Pnext = Pnext + percentStateNoise^2*diag(rho.^2);

end
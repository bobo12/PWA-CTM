function  [rhoNext, Pnext] = EKFaprioriMode(rho, P, percentStateNoise, J, w, d, rhoJ)
%
% rho is an ensemble of vectors of the densities on the road at time j
%

% forecast step using godunov scheme
nExt = length(rho);
s = zeros(nExt-1, 1);
m = zeros(nExt-2, 1);

rhoNext = zeros(nExt, 1);
Pnext = zeros(nExt, nExt);

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
    rhoNext(i) = J(m(i-1),:) * [rho(i-1); rho(i); rho(i+1)] + w(m(i-1));
end

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

Pnext = Pnext + percentStateNoise;

end
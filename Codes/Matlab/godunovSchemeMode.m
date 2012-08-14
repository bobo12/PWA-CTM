function rho_next = godunovSchemeMode(rho, J, w, d, up, dn)
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
    rho_next(i) = J(m(i-1),:) * [rho(i-1); rho(i); rho(i+1)] + w(m(i-1));
end

end
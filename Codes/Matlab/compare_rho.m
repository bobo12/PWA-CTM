function q = compare_rho(rho1,rho2,rho_c,q_rho1,q_rho2,q_rho_c)
if rho_c<=rho2 && rho2<=rho1
    q = q_rho2;
else if rho2<=rho_c && rho_c<=rho1
        q = q_rho_c;
    else if  rho2<=rho1 && rho1<=rho_c
            q = q_rho1;
        else %if rho1<=rho2
                q = min(q_rho1,q_rho2);
            %end
        end
    end
end
end
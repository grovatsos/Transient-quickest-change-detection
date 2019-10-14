function [ fa_Bcusum] = DCUSUMtwoFAB( a0,b0,a1,b1, a2,b2,n,threshold_h,r)

iterations = n;    




threshold = threshold_h;


rho_1_2 = r;


mu_0 = a0;
s_d_0 = b0;

mu_1 = a1;
s_d_1 = b1;

mu_2 = a2;
s_d_2 = b2;

for t=1:1:length(threshold)

    for j =1:1:iterations
        k=1;
        j;
        while 1
            Z(k)=normrnd(mu_0,s_d_0);
            if k == 1 ; 
                Omega_1(k)= log((normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0)));
                Omega_2(k)= log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)));
            else
                Omega_1(k) =   subplus( Omega_1(k-1))+ log((1-r)*(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0))) ;
               Omega_2(k) =  max(Omega_2(k-1)+log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))),Omega_1(k-1) + log(r*(normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)))) ;
            end
            W(k) = subplus(max(Omega_1(k),Omega_2(k)));
            if W(k) > threshold(t)
                 delay(j)=k-1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
            end
            k=k+1;
        end

    end
    BcusumFA(t)=mean(delay)
end
fa_Bcusum=BcusumFA
end

function [FAdcusum]=DCUSMfiveFA(  a0,b0,a1,b1, a2,b2,a3,b3,a4,b4,a5,b5,n,threshold_dcusum,r)
iterations = n;
threshold =threshold_dcusum;
rho_trans = r;
rho_1_2= rho_trans;
rho_2_3= rho_trans;
rho_3_4= rho_trans;
rho_4_5= rho_trans;

mu_0 = a0;
s_d_0 = b0;

mu_1 = a1;
s_d_1 = b1;

mu_2 = a2;
s_d_2 = b2;

mu_3 = a3;
s_d_3 = b3;

mu_4 = a4;
s_d_4 = b4;

mu_5 = a5;
s_d_5 = b5;

delay=0;

for t=1:1:length(threshold)
    for j =1:1:iterations
        %Gamma_1 = geornd(0.2)+1; % Regular value
        % Generating the data
        %Calculating D-CuSum Statistic
        k=1;
        while 1
            Z(k)=normrnd(mu_0,s_d_0);
            if k == 1 ; 
                Omega_1(k)= log((normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0)));
                Omega_2(k)= log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)));
                Omega_3(k)= log((normpdf(Z(k),mu_3,s_d_3))/(normpdf(Z(k),mu_0,s_d_0)));
                Omega_4(k)= log((normpdf(Z(k),mu_4,s_d_4))/(normpdf(Z(k),mu_0,s_d_0)));
                Omega_5(k)= log((normpdf(Z(k),mu_5,s_d_5))/(normpdf(Z(k),mu_0,s_d_0)));
            else
                Omega_1(k) =  subplus(Omega_1(k-1)) + log((normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0))) ;
                Omega_2(k) =  max(Omega_2(k-1),Omega_1(k-1)) + log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))) ;
                Omega_3(k) =  max(Omega_3(k-1),Omega_2(k-1)) + log((normpdf(Z(k),mu_3,s_d_3))/(normpdf(Z(k),mu_0,s_d_0))) ;
                Omega_4(k) =  max(Omega_4(k-1),Omega_3(k-1)) + log((normpdf(Z(k),mu_4,s_d_4))/(normpdf(Z(k),mu_0,s_d_0))) ;
                Omega_5(k) =  max(Omega_5(k-1),Omega_4(k-1)) + log((normpdf(Z(k),mu_5,s_d_5))/(normpdf(Z(k),mu_0,s_d_0))) ;
            end

            
            W_1(k) = max(Omega_1(k),Omega_2(k));
            W_2(k) = max(W_1(k),Omega_3(k));
            W_3(k) = max(W_2(k),Omega_4(k));
            W_4(k) = max(W_3(k),Omega_5(k));
            W(k)=subplus(W_4(k));
            
            if W(k) > threshold(t)
                 delay(j)=k-1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
            end
            k=k+1;
        end
    end
    FAdcusum(t)=mean(delay)
end
FAdcusum
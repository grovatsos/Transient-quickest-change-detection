
function [ add_cusum] = addcusumtwo2( a0,b0,a1,b1, a2,b2,n,threshold_cusum,r,D1)

iterations = n;



threshold =threshold_cusum;


rho_1_2 = r;


mu_0 = a0;
s_d_0 = b0;

mu_1 = a1;
s_d_1 = b1;

mu_2 = a2;
s_d_2 = b2;

for t=1:1:length(threshold)

    clear delay
    clear Gamma_1
    clear Gamma_2
    clear i
    clear j
    clear k
    clear Omera_1
    clear Omega_2
    clear W
    clear Z
    for j =1:1:iterations
        
        



        %Gamma_1 = geornd(0.2)+1; % Regular value
        Gamma_1 = 1;    %Since we are doing Delay here
        Gamma_2 = Gamma_1 + D1;
        horizon=Gamma_2+10000;
        % Generating the data

  %          Z(1:(Gamma_1-1)) = normrnd(mu_0,s_d_0,1,Gamma_1-1);

            Z(Gamma_1:(Gamma_2-1)) = normrnd(mu_1,s_d_1,1,Gamma_2-Gamma_1);

            Z(Gamma_2:horizon) = normrnd(mu_2,s_d_2,1,horizon-Gamma_2+1);


        %Calculating D-CuSum Statistic

        for k = 1:1:horizon
            if k == 1 ; 
                Omega_1(k)= log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)));
                %Omega_2(k)= log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)));
            else
                Omega_1(k) =  subplus(Omega_1(k-1)) + log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))) ;
                %Omega_2(k) =  max(Omega_2(k-1),Omega_1(k-1)) + log((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))) ;
            end

            W(k) = subplus(Omega_1(k));
            if W(k) > threshold(t)
                 delay(j)=k-Gamma_1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
             end
        end


    end
    cusumADD(t)=mean(delay)

end

add_cusum=cusumADD;
end

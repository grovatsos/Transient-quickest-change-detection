function [ add_baysian] = NEWTESTtwotransientperiodsADD( a0,b0,a1,b1, a2,b2,n,threshold_b,r)


iterations = n;

threshold = threshold_b; 


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
        

        Gamma_1 = 1;    %Since we are doing Delay here
        Gamma_2 = Gamma_1 + geornd(rho_1_2) +1;
        horizon = Gamma_2+10000;
        % Generating the data

  %      for i =1:1:(Gamma_1-1)
            Z(1:Gamma_1-1) = normrnd(mu_0,s_d_0,1,Gamma_1-1);
 %       end
 %       for i =Gamma_1:1:(Gamma_2-1)
            Z(Gamma_1:(Gamma_2-1)) = normrnd(mu_1,s_d_1,1,(Gamma_2)-Gamma_1);
 %       end
 %       for i =Gamma_2:1:horizon
            Z(Gamma_2:horizon) = normrnd(mu_2,s_d_2,1,horizon-Gamma_2+1);
 %       end

        %Calculating D-CuSum Statistic
        for k = 1:1:horizon
            if k == 1 ; 
                q_1(k)=(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0));
                q_2(k)=0;
            else
                q_1(k)=(1+q_1(k-1)*(1-rho_1_2))*((normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0)));
                q_2(k)=(q_1(k-1)*rho_1_2 +q_2(k-1))*((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))); 
            end
            

            
            W(k) = log(q_1(k)+q_2(k));
            if W(k) > threshold(t)
                 delay(j)=k-Gamma_1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
             end
        end


    end
    baysianADD(t)=mean(delay)

end

add_baysian=baysianADD
end

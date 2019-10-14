function [PFA, ADD] = PFAADDBayes( a0,b0,a1,b1, a2,b2,n,threshold_b,r,r1)
iterations = n;

threshold = threshold_b; 
% horizon = 100000;

rho_1_2 = r;
rho=r1;

mu_0 = a0;
s_d_0 = b0;

mu_1 = a1;
s_d_1 = b1;

mu_2 = a2;
s_d_2 = b2;
delay=0;
ADD=0;
for t=1:1:length(threshold)
    FA=0;
    t
    for j =1:1:iterations
%         if mod(j,100) ==0
%             j
%         end
        %Gamma_1 = geornd(0.2)+1; % Regular value
        Gamma_1 = geornd(rho)+1;    %Since we are doing Delay here
        Gamma_2 = Gamma_1 + geornd(rho_1_2) +1;
        horizon=Gamma_2+10000;
        % Generating the data

        Z(1:Gamma_1-1) = normrnd(mu_0,s_d_0,1,Gamma_1-1);
 
        Z(Gamma_1:(Gamma_2-1)) = normrnd(mu_1,s_d_1,1,(Gamma_2)-Gamma_1);

        Z(Gamma_2:horizon) = normrnd(mu_2,s_d_2,1,horizon-Gamma_2+1);

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
%                 delay(j)=k-Gamma_1;
                 if k<Gamma_1
                     FA=FA+1;
                     break;
                 else
                     delay(j)=k-Gamma_1;
                     break
                 end
             end
        end


    end
    baysianADD(t)=sum(delay)/(iterations-FA)
    PFA(t)=FA/iterations
end
PFA
ADD=baysianADD
end
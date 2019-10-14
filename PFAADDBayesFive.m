function [ PFA,ADD ] = PFAADDBayesFive(  a0,b0,a1,b1, a2,b2,a3,b3,a4,b4,a5,b5,n,threshold_b,r,r1)



iterations = n;




threshold = threshold_b;


rho = r1; 

rho_1_2= r;
rho_2_3= r;
rho_3_4= r;
rho_4_5= r;

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
t
    FA=0;
    for j =1:1:iterations

%         if mod(j,100) ==0
%             j
%         end
        
        %Gamma_1 = geornd(0.2)+1; % Regular value
        Gamma_1 = geornd(rho)+1;    %Since we are doing Delay here
        Gamma_2 = Gamma_1 + geornd(rho_1_2) +1;
        Gamma_3 = Gamma_2 + geornd(rho_2_3) +1;
        Gamma_4 = Gamma_3 + geornd(rho_3_4) +1;
        Gamma_5 = Gamma_4 + geornd(rho_4_5) +1;
        horizon=Gamma_5+10000;

        % Generating the data


            Z(1:(Gamma_1-1)) = normrnd(mu_0,s_d_0,1,Gamma_1-1);

            Z(Gamma_1:(Gamma_2-1)) = normrnd(mu_1,s_d_1,1,Gamma_2-Gamma_1);


            Z(Gamma_2:(Gamma_3-1)) = normrnd(mu_2,s_d_2,1,Gamma_3-Gamma_2);

            Z(Gamma_3:Gamma_4-1) = normrnd(mu_3,s_d_3,1,Gamma_4-Gamma_3);

            Z(Gamma_4:Gamma_5-1) = normrnd(mu_4,s_d_4,1,Gamma_5-Gamma_4);

            Z(Gamma_5:horizon) = normrnd(mu_5,s_d_5,1,horizon-Gamma_5+1);

        %Calculating New Test Statistic

        for k = 1:1:horizon
            if k == 1 ; 
                q_1(k)=(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0));
                q_2(k)=0;
                q_3(k)=0;
                q_4(k)=0;
                q_5(k)=0;
            else
                q_1(k)=(1+q_1(k-1)*(1-rho_1_2))*((normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0)));
                q_2(k)=(q_1(k-1)*rho_1_2 +q_2(k-1)*(1-rho_2_3))*((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))); 
                q_3(k)=(q_2(k-1)*rho_2_3 +q_3(k-1)*(1-rho_3_4))*((normpdf(Z(k),mu_3,s_d_3))/(normpdf(Z(k),mu_0,s_d_0))); 
                q_4(k)=(q_3(k-1)*rho_3_4 +q_4(k-1)*(1-rho_4_5))*((normpdf(Z(k),mu_4,s_d_4))/(normpdf(Z(k),mu_0,s_d_0))); 
                q_5(k)=(q_4(k-1)*rho_4_5 +q_5(k-1))*((normpdf(Z(k),mu_5,s_d_5))/(normpdf(Z(k),mu_0,s_d_0))); 
            end
            

            
            W(k) = log(q_1(k)+q_2(k)+q_3(k)+q_4(k)+q_5(k));
            if W(k) > threshold(t)
                if k<Gamma_1
                    FA=FA+1;
                    break
                else
                 delay(j)=k-Gamma_1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
                end
%             else
%                 delay=100000000000000000;
             end
        end

    end
    pfa(t)=FA/iterations;
    add(t)=sum(delay)/(iterations-FA);
end
PFA=pfa
ADD=add




end


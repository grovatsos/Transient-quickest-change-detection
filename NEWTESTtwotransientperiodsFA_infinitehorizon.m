clear all
clc

iterations = 1000;




threshold = [1:9];
%horizon = 10000;

rho = 0.2; 
rho_1_2 = 0.1;


mu_0 = 0;
s_d_0 = 2;

mu_1 = 1;
s_d_1 = 2;

mu_2 = 2;
s_d_2 = 2;

for t=1:1:length(threshold)
threshold(t)
    for j =1:1:iterations

%         if mod(j,100) ==0
%             j
%         end
        
%    Z=normrnd(mu_0,s_d_0,1,horizon);

        %Calculating New Test Statistic
        k=1;
        
        while 1
            Z(k)=normrnd(mu_0,s_d_0);
            if k == 1 ; 
                q_1(k)=(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0));
                q_2(k)=0;
            else
                q_1(k)=(1+q_1(k-1)*(1-rho_1_2))*((normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0)));
                q_2(k)=(q_1(k-1)*rho_1_2 +q_2(k-1))*((normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))); 
            end
            

            
            W(k) = log(q_1(k)+q_2(k));
            if W(k) > threshold(t)
                 delay(j)=k-1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
            end
            k=k+1;
        end

    end
    FA(t)=mean(delay)
end
FA

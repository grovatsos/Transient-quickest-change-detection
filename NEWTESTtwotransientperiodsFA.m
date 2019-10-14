function [ fa_baysian] = NEWTESTtwotransientperiodsFA( a0,b0,a1,b1, a2,b2,n,threshold_b,r)

iterations = n;




threshold =threshold_b;



rho_1_2 = r;


mu_0 = a0;
s_d_0 = b0;

mu_1 = a1;
s_d_1 = b1;

mu_2 = a2;
s_d_2 = b2;

for t=1:1:length(threshold)
    clear delay
 
    for j =1:1:iterations

%         if mod(j,100) ==0
%             j
%         end
        


        %Calculating New Test Statistic
        k=1;   
        while 1
                Z(k)=normrnd(mu_0,s_d_0,1);
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
    baysianFA(t)=mean(delay)
end
fa_baysian=baysianFA
end
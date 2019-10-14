
function[FAbayesian]=NEWTESTfivetransientperiodsFA(  a0,b0,a1,b1, a2,b2,a3,b3,a4,b4,a5,b5,n,threshold_b,r)
iterations = n;
threshold = threshold_b;
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
       %Calculating New Test Statistic
        k=1;
        while 1
            Z(k)=normrnd(mu_0,s_d_0);
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
                 delay(j)=k-1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
            %else
             %   delay(j)=100000000000000000;
            end
             k=k+1;
        end

    end
    FAbayesian(t)=mean(delay)
end
FAbayesian
function [ fa_h] = HFA( a0,b0,a1,b1, a2,b2,n,threshold_h,r)

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

%         if mod(j,100) ==0
%            j
%         end

        k=1;
        
        while 1
            Z(k)=normrnd(mu_0,s_d_0);
             if k == 1 ; 
                h_1(k)= 0;
                h_2(k)= (1-rho_1_2)*(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0));
            else
                h_1(k) =  h_1(k-1)*(normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))+rho_1_2/(1-rho_1_2) *  (normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)) * h_2(k-1);
                h_2(k) =  (1-rho_1_2)*(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0))*(1+h_2(k-1));
            end
            W(k) = log(h_1(k)+h_2(k));
            if W(k) > threshold(t)
                 delay(j)=k-1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
            end
            k=k+1;
        end

    end
    hfa(t)=mean(delay)
end
fa_h=hfa
end

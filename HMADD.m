
function [ add_hm] = HMADD( a0,b0,a1,b1, a2,b2,n,threshold_h,r)

iterations = n;



threshold =threshold_h;


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
        Gamma_2 = Gamma_1 + geornd(rho_1_2) +1;
        horizon=Gamma_2+10000;
        % Generating the data

  %          Z(1:(Gamma_1-1)) = normrnd(mu_0,s_d_0,1,Gamma_1-1);

            Z(Gamma_1:(Gamma_2-1)) = normrnd(mu_1,s_d_1,1,Gamma_2-Gamma_1);

            Z(Gamma_2:horizon) = normrnd(mu_2,s_d_2,1,horizon-Gamma_2+1);


        %Calculating D-CuSum Statistic

        for k = 1:1:horizon
             if k == 1 ; 
                h_1(k,1)= 0;
                h_2(k,1)= (1-rho_1_2)*(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0));
                h= h_1(k,:)+h_2(k,:);
             else
                 h_1(k,1:(k-1)) =  h_1(k-1,1:(k-1))*(normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0))+rho_1_2/(1-rho_1_2) *  (normpdf(Z(k),mu_2,s_d_2))/(normpdf(Z(k),mu_0,s_d_0)) * h_2(k-1,1:(k-1));
                 h_2(k,1:(k-1)) =  (1-rho_1_2)*(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0))*(h_2(k-1,1:(k-1)));

                 h_1(k,k)=0;
                 h_2(k,k)=(1-rho_1_2)*(normpdf(Z(k),mu_1,s_d_1))/(normpdf(Z(k),mu_0,s_d_0));
                 h=h_1(k,:)+h_2(k,:);
            end
            W(k) = log(max(h));
            if W(k) > threshold(t)
                 delay(j)=k-Gamma_1;
                 %fprintf('Crossed a threshold of %d. at time instant %d. with a delay of  %d.',threshold, k,delay(j))
                 break
             end
        end


    end
    hmADD(t)=mean(delay)

end

add_hm=hmADD;
end

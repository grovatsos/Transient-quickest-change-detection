clear all
clc
tic
rho_trans = 0.02;


mu_0 = 0;
s_d_0 = 1;

mu_1 = 3;
s_d_1 = 1;

mu_2 = 3;
s_d_2 = 1;

mu_3 = 5;
s_d_3 =1;

mu_4 =3;
s_d_4 = 1;

mu_5 = 0;
s_d_5 =0.9;

threshold_b=[0.1:1:5.1];%[0.1:0.9:7,7.1,7.2];
threshold_cusum=[0.01,0.02,0.1:0.7:5];%[0.1,0.5,0.7,0.9:1.1:6];


rho_1_2= rho_trans;
rho_2_3= rho_trans;
rho_3_4= rho_trans;
rho_4_5= rho_trans;
n=100%;%number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fa_dcusum=DCUSMfiveFA(  mu_0,s_d_0,mu_1,s_d_1, mu_2,s_d_2,mu_3,s_d_3,mu_4,s_d_4,mu_5,s_d_5,n,threshold_cusum,rho_trans);
add_dcusum=DCUSMfiveADD(  mu_0,s_d_0,mu_1,s_d_1, mu_2,s_d_2,mu_3,s_d_3,mu_4,s_d_4,mu_5,s_d_5,n,threshold_cusum,rho_trans);
fa_baysian=NEWTESTfivetransientperiodsFA(  mu_0,s_d_0,mu_1,s_d_1, mu_2,s_d_2,mu_3,s_d_3,mu_4,s_d_4,mu_5,s_d_5,n,threshold_b,rho_trans);
add_baysian= NEWTESTfivetransientperiodsADD(  mu_0,s_d_0,mu_1,s_d_1, mu_2,s_d_2,mu_3,s_d_3,mu_4,s_d_4,mu_5,s_d_5,n,threshold_b,rho_trans);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
loglog(fa_dcusum,add_dcusum,'-.rs','LineWidth',2,'MarkerSize',8);
hold;
loglog(fa_baysian,add_baysian,'--g*','LineWidth',2,'MarkerSize',8);
legend('D-CuSum','Bayesian');
xlabel('Mean Time to False Alarm');
ylabel('Detection Delay');

%%%%%%%%%


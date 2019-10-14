%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
tic
rho_1_2 = 0.25;


mu_0 = 0;
s_d_0 = 1;

mu_1 = 1;
s_d_1 =1;

mu_2 =2;
s_d_2 = 1;
n=100; %iterations
threshold_h=[0.01:0.5:5];
threshold_hm=[0.01:0.5:5];
threshold_b=[4:2:12];%[0.01:0.5:5];%[0:2:20];%[0.1:0.9:7,7.1,7.2];
threshold_cusum=[4:2:12];%[0.1:0.5:6];%[0.1,0.5,0.7,0.9:1.1:6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[fa_baysian]=DCUSUMtwoFAB(mu_0,s_d_0,mu_1 ,s_d_1,mu_2,s_d_2,n,threshold_b,rho_1_2);
[add_baysian]=DCUSMtwoADDB(mu_0,s_d_0,mu_1 ,s_d_1,mu_2,s_d_2,n,threshold_b,rho_1_2);
[add_dcusum]=DCUSMtwotransientperiodsADD(mu_0,s_d_0,mu_1 ,s_d_1,mu_2,s_d_2,n,threshold_cusum,rho_1_2);
[fa_dcusum]=DCUSUMtwotransientperiodsFA(mu_0,s_d_0,mu_1 ,s_d_1,mu_2,s_d_2,n,threshold_cusum,rho_1_2);

toc
figure
semilogx(fa_dcusum,add_dcusum,'-.rs','LineWidth',2,'MarkerSize',8);
hold;
semilogx(fa_baysian,add_baysian,'--g*','LineWidth',2,'MarkerSize',8);
% loglog(hfa,hadd,'--b*','LineWidth',2,'MarkerSize',8);
% loglog(hmfa,hmadd,'--ko','LineWidth',2,'MarkerSize',8);
legend('D-CuSum','Bayesian');
xlabel('Mean Time to False Alarm');
ylabel('WADD');


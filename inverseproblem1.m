% Covid model inverse problem, NETHERLANDS thetas

clear all; close all;

load('data_Ne.mat');
load('T_Ne.mat');
thetas_real = [0.32 -0.38 0.4 -0.37 0.42 -0.4 0.34 -0.25 0.13 -0.0098];
time = T;
data = Inew_real;

%thetas_IC = [1, -1, 0.5, -0.5, 3, -1, 1.75, -1.75,2.1,-2.1];
thetas_IC = zeros(1,10);
thetas_IC(1) = unifrnd(.1,1);


[thetas_fit,resnorm] = lsqcurvefit(@forwardmodel1,thetas_IC,time,data,[],[],...
                             optimset('Disp','iter','TolX',10^(-15),'TolFun',10^(-15)));

error = abs(thetas_fit - thetas_real);
resnorm

figure
plot(error,'*', 'LineWidth',2)
hold on
plot(zeros(1,0),':', 'LineWidth',2)
xlabel(['Leg coef'], 'fontsize',24)
ylabel('error Leg coef', 'fontsize',24)
hold off

figure
plot(time, betacomp(time,10,thetas_fit,0,60),time, betacomp(time,10,thetas_real,0,60),time, betacomp(time,10,thetas_IC,0,60), 'LineWidth',2)
xlabel('time', 'fontsize',24)
ylabel('beta', 'fontsize',24)
lgd = legend('beta_{fit}','beta_{real}', 'beta_{IC}', 'fontsize',24);


figure
plot(time,forwardmodel1(thetas_fit,time), time,data,'*', 'LineWidth',2)
hold on
plot(time,forwardmodel1(thetas_real,time),'-.', 'LineWidth',2)
hold off
  
n_boots =100;
thetas_i = zeros(n_boots,10);
for i =1:n_boots
    i
    data_i = poissrnd(data);
    [thetas_fit,resnorm] = lsqcurvefit(@forwardmodel1,thetas_IC,time,data_i,[],[],...
                             optimset('TolX',10^(-15),'TolFun',10^(-15)));
   thetas_i(i,:)= thetas_fit';
end

figure
histogram(thetas_i(:,1))

function Inew = forwardmodel1(thetas,time)
% params and IC for ODE solver
N = 3*10^8;
vac_perc = 0.5;
IFR_S = 0.012;
IFR_V = IFR_S/10;
end_day=60;

par(1) = N;
par(2) = (0.67-0.56)/end_day; % p
par(3) = 0.75; % alpha 
par(4) = (1-IFR_S) * 1/5; % gammaSR 
par(5) = IFR_S * 1/5; % gammaSD  
par(6) = (1- IFR_V)*1/3; % gammaVR 
par(7) = IFR_V*1/3; % gammaVD 

V0 = N*vac_perc; %1.8*10^8;
I0 = 1*10^5;
Is0 = I0*(1-vac_perc);
Iv0 = I0*vac_perc;
D0 = 6*10^5;
S0 = N-V0-Is0-Iv0-D0;
% ODE solver
[~, y] = ode45(@(t,y) RHS(t,y,thetas,par),time,[S0 V0 Is0 Iv0 D0],[thetas,par]);
% get Inew
beta = betacomp(time,10,thetas,0,60); 
Inew = beta .* y(:,1).* ((y(:,3) +y(:,4)) ./ (N - y(:,5))) + beta.* (1-par(3)) .* y(:,2).* ((y(:,3) +y(:,4)) ./ (N - y(:,5)));
end



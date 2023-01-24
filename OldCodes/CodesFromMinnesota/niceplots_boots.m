% Nice plots, run after inverseproblem1.m
close all;

colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]; ...
        [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]; [0.6350, 0.0780, 0.1840]; [0.8350, 0.6780, 0.1840]];
load('thetas_i.mat')

close all;
t = tiledlayout(2,5);
for i = 1:10
ax1 = nexttile;
histogram(thetas_i(:,i))
xlabel('\beta' +string(i))
ax1.FontSize = 15;
ax1.LineWidth = 1;
end
t.Padding = 'none';
t.TileSpacing = 'none';
saveas(gcf,'betas_uncertainty.png')



close all;
plot(time, betacomp(time,10,thetas_fit,0,60),...
    time,betacomp(time,10,thetas_real,0,60),time, betacomp(time,10,thetas_IC,0,60), 'LineWidth',2)
xlabel('time', 'fontsize',24)
ylabel('beta', 'fontsize',24)
legend('beta_{fit}','beta_{true}', 'beta_{IC}', 'fontsize',24);
set(gca,'Fontsize',20, 'LineWidth',1);
saveas(gcf,'betas_lsqfit.png')


close all;
plot(time,forwardmodel1(thetas_fit,time),time,forwardmodel1(thetas_real,time),'-.', time,data,'b*', 'LineWidth',2)
xlabel('time', 'fontsize',24)
ylabel('cases', 'fontsize',24)
legend('Fit I_{new}', 'True I_{new}', 'Synth data', 'fontsize',24);
set(gca,'Fontsize',20, 'LineWidth',1);
saveas(gcf,'Inew_lsqfit.png')


close all;
for i =1:100
plot(time,forwardmodel1(thetas_i(i,:),time), 'LineWidth',2)
hold on
end
plot(time,forwardmodel1(thetas_real,time),'-.', time,data,'b*', 'LineWidth',2)
xlabel('time', 'fontsize',24)
ylabel('cases', 'fontsize',24)
set(gca,'Fontsize',20, 'LineWidth',1);
saveas(gcf,'Inew_uncertainty.png')


function Inew = forwardmodel1(thetas,time)
% params and IC for ODE solver
[N, vac_perc, IFR_S, IFR_V, end_day, par] = parameters_setup();
[S0, V0, Is0, Iv0, D0] = initial_condition_setup(N, vac_perc);
% ODE solver
[~, y] = ode45(@(t,y) RHS(t,y,thetas,par),time,[S0 V0 Is0 Iv0 D0],[thetas,par]);
% get Inew
beta = betacomp(time,10,thetas,0,60); 
Inew = beta .* y(:,1).* ((y(:,3) +y(:,4)) ./ (N - y(:,5))) + beta.* (1-par(3)) .* y(:,2).* ((y(:,3) +y(:,4)) ./ (N - y(:,5)));
end

function [N, vac_perc, IFR_S, IFR_V, end_day,par] = parameters_setup()
% params 
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
end

function [S0, V0, Is0, Iv0, D0] = initial_condition_setup(N, vac_perc)
V0 = N*vac_perc; %1.8*10^8;
I0 = 1*10^5;
Is0 = I0*(1-vac_perc);
Iv0 = I0*vac_perc;
D0 = 6*10^5;
S0 = N-V0-Is0-Iv0-D0;
end
% Covid model forward simulation, NETHERLANDS thetas

clear all; close all;
% Fix the seed to have reproducible data:
seed = 10;
rng(seed);

% initial conditions
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

% solve differential equations for vector y from day 0 to day end_day
% and initial conditions as specified


thetas = [0.32 -0.38 0.4 -0.37 0.42 -0.4 0.34 -0.25 0.13 -0.0098];
%thetas = [0.43 -0.73 0.92 -1 1.1 -0.97 0.92 -0.59 0.31 -0.08];
[T, y] = ode45(@(t,y) RHS(t,y,thetas,par),[0 end_day],[S0 V0 Is0 Iv0 D0],[thetas,par]);


Itotal = y(:,3) + y(:,4);



beta = betacomp(T,10,thetas,0,60); 
Inew = beta .* y(:,1).* ((y(:,3) +y(:,4)) ./ (N - y(:,5))) + beta.* (1-par(3)) .* y(:,2).* ((y(:,3) +y(:,4)) ./ (N - y(:,5)));
Inew_real = Inew + randn(size(Inew))*10^5;
for i =1:length(Inew)
    if Inew_real(i)<0
        Inew_real(i) = Inew(i) + unifrnd(-Inew(i),Inew(i));
    end
end
figure
plot(T,Inew,T,Inew_real,'*', 'LineWidth',2)
title('Daily infections')

D_daily = y(2:end,5) - y(1:end-1,5);
D_daily_noise = D_daily +randn(size(D_daily))*10^2;
for i =1:length(D_daily)
    if D_daily_noise(i)<0
        D_daily_noise(i) = D_daily(i) + unifrnd(-D_daily(i),D_daily(i));
    end
end
figure
plot(T(1:end-1),D_daily, 'LineWidth',2)
hold on
plot(T(1:end-1),D_daily_noise,'*', 'LineWidth',2)
D_real = y(1,5) + [1 ;cumsum(D_daily_noise)];
figure
plot(T,y(:,5))
hold on
plot(T,D_real,'*')


size(D_real)
size(Inew_real)

save('data_Ne.mat','Inew_real','T', 'D_real')
%save('T_Ne.mat','T')

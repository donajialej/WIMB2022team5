% SIR model from class Jan 4th

clear all;

% initial conditions
N = 3*10^8;
V0 = 1.8*10^8;
Is0= 5*10^4;
Iv0 = 5*10^4;
D0 = 6*10^5;
S0 = N-V0-Is0-Iv0-D0;

beta = 1;

par(1) = N;
par(2) = 0.01; % p
par(3) = 0.75; % alpha 
par(4) = (1-0.018) * 1/5; % gammaSR 
par(5) = 0.018 * 1/5; % gammaSD  
par(6) = (1- 0.006)*1/3; % gammaVR 
par(7) = 0.006*1/3; % gammaVD 

% solve differential equations for vector y=[S;I] from day 0 to day end_day
% and initial conditions as specified

end_day=60;
time = linspace(0,60);
thetas = [0.32 -0.38 0.4 -0.37 0.42 -0.4 0.34 -0.25 0.13 -0.0098];

[T, y] = ode45(@(t,y) RHS(t,y,par,thetas),[0 end_day],[S0 V0 Is0 Iv0 D0],[par,thetas]);

plot(T,y(:,3:5), 'LineWidth',2)
lgd = legend('I_s','I_v','D');

Itotal = y(:,3) + y(:,4);
Inew = Itotal(2:end)-Itotal(1:end-1);
figure
plot(T,y(:,3) + y(:,4), 'LineWidth',2)
title('Total infections')

figure
plot(T(2:end),y(2:end,5) - y(1:end-1,5), 'LineWidth',2)
title('Daily deaths')

figure
plot(T(2:end),Inew, 'LineWidth',2)
title('Daily infections')


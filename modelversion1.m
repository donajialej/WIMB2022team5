% SIR model from class Jan 4th

clear all; close all;

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

% solve differential equations for vector y=[S;I] from day 0 to day end_day
% and initial conditions as specified


time = linspace(0,60);
thetas = [0.32 -0.38 0.4 -0.37 0.42 -0.4 0.34 -0.25 0.13 -0.0098];

[T, y] = ode45(@(t,y) RHS(t,y,par,thetas),[0 end_day],[S0 V0 Is0 Iv0 D0],[par,thetas]);

beta = betacomp(T,10,thetas,0,60); 

plot(T,y(:,3:5), 'LineWidth',2)
lgd = legend('I_s','I_v','D');

Itotal = y(:,3) + y(:,4);
Inew = beta .* y(:,1).* ((y(:,3) +y(:,4)) ./ (N - y(:,5))) + beta.* (1-par(3)) .* y(:,2).* ((y(:,3) +y(:,4)) ./ (N - y(:,5)));
Inew_real = Inew + randn(size(Inew))*10^4;
for i =1:length(Inew)
    if Inew_real(i)<0
        Inew_real(i) = Inew(i) + unifrnd(-Inew(i),Inew(i));
    end
end


figure
plot(T,y(:,3) + y(:,4), 'LineWidth',2)
title('Total infections')

figure
plot(T(2:end),y(2:end,5) - y(1:end-1,5), 'LineWidth',2)
title('Daily deaths')

figure
plot(T,Inew,T,Inew_real,'*', 'LineWidth',2)
title('Daily infections')


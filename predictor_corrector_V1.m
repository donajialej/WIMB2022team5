% testing my fourier expansions
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
par(8) = end_day;

V0 = N*vac_perc; %1.8*10^8;
I0 = 1*10^5;
Is0 = I0*(1-vac_perc);
Iv0 = I0*vac_perc;
D0 = 6*10^5;
S0 = N-V0-Is0-Iv0-D0;
thetas = [0.32 -0.38 0.4 -0.37 0.42 -0.4 0.34 -0.25 0.13 -0.0098];


n_u = [2,2,2,2,2];

u = ones(2,10);

dy = forward_expansion(10,u,thetas,par, n_u);
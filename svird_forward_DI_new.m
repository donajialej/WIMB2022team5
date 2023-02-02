function svird_forward_DI_new

close all 
clear all
clc  
format long 
warning('off','all')

global N alpha p gammasd gammasr gammavr gammavd delta1 delta2 m theta a b

a = 1;
b = 140;
m = 10;
tdata = a:1:b;

alpha = .8;
% .8
gammavd = 0.005/12.7/18.5;   
gammasd = 0.005/18.5;
gammasr = (1-0.005)/10;
gammavr = (1-0.005/12.7)/10;
p = 0.00086;

delta1 = 1/90;
% 90
delta2 = 0;
% 0

% Total population of the region

N = 10799566; % GA population

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta(1)=0.26; theta(2)=-0.13; theta(3)=0.093; theta(4)=-0.036; theta(5)=0.05;
theta(6)=-0.12; theta(7)=0.14; theta(8)=-0.098; theta(9)=0.031; theta(10)=0.023;

% Y0(1) = (N - 3942002 - 3580 - (1139396-21486))/N;
Y0(1) = N - 3942002 - 3580;
Y0(2) = 3942002;

Y0(3) = 3580 - 731;
Y0(4) = 731;
Y0(5) = 0;
Y0(6) = 0;

[~,Y] = ode23s(@(t,y) svird(t,y),tdata,Y0);

Idata_clean = beta(tdata)'.*Y(:,1).*(Y(:,3)+Y(:,4))./(N - Y(:,5)) +...
    (1-alpha)*beta(tdata)'.*Y(:,2).*(Y(:,3)+Y(:,4))./(N - Y(:,5));

Ddata_clean = gammasd.*Y(:,3) + gammavd.*Y(:,4);

%%%%%%%%%%%%%%%%%% GENERATE ACTUAL INCIDENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Idata = Idata_clean;
    Ddata = Ddata_clean;
    for i = 1:length(tdata)
    Idata(i) = Idata_clean(i) + normrnd(0,20000);
    Ddata(i) = Ddata_clean(i) + normrnd(0,30);
    if Idata(i) < 0
       Idata(i) = Idata_clean(i) + unifrnd(-Idata_clean(i),Idata_clean(i));
    end 
    if Ddata(i) < 0
       Ddata(i) = Ddata_clean(i) + unifrnd(-Ddata_clean(i),Ddata_clean(i));
    end 
    end 
    
    D = cumsum(Ddata) + Y0(6);
    D_clean = cumsum(Ddata_clean) + Y0(6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure; 
axes2 = axes('Parent',figure1,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all');

line1 = plot(tdata,Idata_clean,'r');
set(line1,'Linewidth',2)
hold on

line2 = plot(tdata,Idata,'*k');
set(line2,'LineWidth',2)

legend([line1(1) line2(1)],{'True Incidence Curve', 'Incidence Curve with Noise Added'},'FontSize',12,'Location','best');

xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14}Case incidence')

figure(figure1)

figure4 = figure; 
axes2 = axes('Parent',figure4,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all');

line11 = plot(tdata,Ddata_clean,'r');
set(line11,'Linewidth',2)
hold on

line22 = plot(tdata,Ddata,'*k');
set(line22,'LineWidth',2)

legend([line11(1) line22(1)],{'True Daily Deaths', 'Daily Deaths with Noise Added'},'FontSize',12,'Location','best');

xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14}Daily Deaths')

figure(figure4)

figure5 = figure; 
axes2 = axes('Parent',figure5,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all');

line111 = plot(tdata,D_clean,'r');
set(line111,'Linewidth',2)
hold on

line222 = plot(tdata,D,'*k');
set(line222,'LineWidth',2)

legend([line111(1) line222(1)],{'True Cumulative Deaths', 'Cumulative Deaths with Noise Added'},'FontSize',12,'Location','best');

xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14}Cumulative Deaths')

figure(figure5)

figure2 = figure; 
axes2 = axes('Parent',figure2,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all');

line3 = plot(tdata,beta(tdata),'b');
set(line3,'Linewidth',2)

legend(line3(1),{'Transmision Rate'},'FontSize',12,'Location','best');

xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14}Transmission rate')

figure(figure2)

figure3 = figure;
subplot(3,2,1)
line77 = plot(tdata,Y(:,1), '-g');
legend(line77(1),{'Reconstructed S(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} S(t)')
% hold on 
subplot(3,2,2)
line78 = plot(tdata,Y(:,2), 'Color', [0.6350, 0.0780, 0.1840]);
legend(line78(1),{'Reconstructed V(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} V(t)')
% hold on 
subplot(3,2,3)
line79 = plot(tdata,Y(:,3), '-m');
legend(line79(1),{'Reconstructed I_s(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} I_s(t)')
% hold on 
subplot(3,2,4)
line76 = plot(tdata,Y(:,4), '-c');
legend(line76(1),{'Reconstructed I_v(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} I_v(t)')
% hold on 
subplot(3,2,5)
line767 = plot(tdata,Y(:,5), '-r');
legend(line767(1),{'Reconstructed R(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} R(t)')
figure(figure3)

SynthDI1R = [tdata' Idata Ddata D];
writematrix(SynthDI1R,'SynthDI1R.txt','Delimiter','tab')
%writematrix(SynthDI1R)

end 

% FUNCTION DEFINITIONS

function dydt = svird(t,y)
global N alpha p gammasd gammasr gammavr gammavd delta1 delta2
dydt = zeros(6,1);
dydt(1) = -beta(t).*y(1).*(y(3)+y(4))./(N - y(6)) - p.*y(1) + delta1*y(5) + delta2*y(2);
dydt(2) = p.*y(1) - (1-alpha)*beta(t).*y(2).*(y(3)+y(4))./(N - y(6)) - delta2*y(2);
dydt(3) = beta(t).*y(1).*(y(3)+y(4))./(N - y(6)) - (gammasr + gammasd)*y(3);
dydt(4) = (1-alpha).*beta(t).*y(2).*(y(3)+y(4))./(N - y(6)) - (gammavr + gammavd)*y(4);
dydt(5) = gammasr.*y(3) + gammavr.*y(4) - delta1*y(5);
dydt(6) = gammasd.*y(3) + gammavd.*y(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yy = beta(t) 
global m theta
       yy = 0;
   for j = 1:m
       yy = yy + theta(j).*leg(j-1,t);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = leg(j,t)
global a b
x = (2.*t - a - b)./(b - a);

if j == 0
   P1 = 1; P = P1; 
elseif j == 1 
   P2 = x; P = P2; 
else   
   P1 = 1; P2 = x; 
       for k = 2:j
           P3 = ((2*(k-1)+1).*x.*P2 - (k-1).*P1)./k;
           P1 = P2; P2 = P3;   
       end
           P = P3;
end
end



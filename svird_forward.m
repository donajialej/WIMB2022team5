function svird_forward

close all 
clear all
clc  
format long 
warning('off','all')

global N alpha p gammasd gammasr gammavr gammavd m theta a b

a = 0;
b = 90;
m = 10;

tdata = a:1:b;
alpha = .75;
gammasd = .012*1/5;
gammavd = .0012*1/4;
gammasr = (1-.012)*1/5;
gammavr = (1-.0012)*1/4;
p = (.67-.56)/(b-a);

N = 10620000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta(1)=0.26; theta(2)=-0.13; theta(3)=0.093; theta(4)=-0.036; theta(5)=0.05;
theta(6)=-0.12; theta(7)=0.14; theta(8)=-0.098; theta(9)=0.031; theta(10)=0.023;

vac_perc = .35;
inf_perc = .0003;
Y0(1) = (1 - vac_perc)*(1-inf_perc)*N;
Y0(2) = vac_perc*(1-inf_perc)*N;
I0 = inf_perc*N;

Y0(3) = I0*(1-vac_perc);
Y0(4) = I0*vac_perc;
Y0(5) = 0.012*Y0(3)+0.0012*Y0(4);

[~,Y] = ode45(@(t,y) svird(t,y),tdata,Y0);

Idata_clean = beta(tdata)'.*Y(:,1).*(Y(:,3)+Y(:,4))./(N - Y(:,5)) +...
    (1-alpha)*beta(tdata)'.*Y(:,2).*(Y(:,3)+Y(:,4))./(N - Y(:,5));

size(Idata_clean)
%%%%%%%%%%%%%%%%%%GENERATE ACTUAL INCIDENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Idata = Idata_clean;
    for i = 1:length(tdata)
    Idata(i) = Idata_clean(i) + normrnd(0,500);
    if Idata(i) < 0
       Idata(i) = Idata_clean(i) + unifrnd(-Idata_clean(i),Idata_clean(i));
    end 
    end 
    
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

figure2 = figure; 
axes2 = axes('Parent',figure2,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all');

line3 = plot(tdata,beta(tdata),'b');
set(line3,'Linewidth',2)

legend([line3(1)],{'Transmision Rate'},'FontSize',12,'Location','best');

xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14}Case incidence')

figure(figure2)

figure3 = figure;
subplot(2,2,1)
line77 = plot(tdata,Y(:,1), '-g');
legend([line77(1)],{'Reconstructed S(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} S(t)')
% hold on 
subplot(2,2,2)
line78 = plot(tdata,Y(:,2), 'Color', [0.6350, 0.0780, 0.1840]);
legend([line78(1)],{'Reconstructed V(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} V(t)')
% hold on 
subplot(2,2,3)
line79 = plot(tdata,Y(:,3), '-m');
legend([line79(1)],{'Reconstructed I_s(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} I_s(t)')
% hold on 
subplot(2,2,4)
line76 = plot(tdata,Y(:,4), '-c');
legend([line76(1)],{'Reconstructed I_v(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} I_v(t)')
figure(figure3)

end 

% FUNCTION DEFINITIONS

function dydt = svird(t,y)
global N alpha p gammasd gammasr gammavr gammavd 
dydt = zeros(5,1);
dydt(1) = beta(t).*y(1).*(y(3)+y(4))./(N - y(5)) - p.*y(1);
dydt(2) = p.*y(1) - (1-alpha)*beta(t).*y(2).*(y(3)+y(4))./(N - y(5));
dydt(3) = beta(t).*y(1).*(y(3)+y(4))./(N - y(5)) - (gammasr + gammasd)*y(3);
dydt(4) = (1-alpha).*beta(t).*y(2).*(y(3)+y(4))./(N - y(5)) - (gammavr + gammavd)*y(4);
dydt(5) = gammasd.*y(3) + gammavd.*y(4);
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



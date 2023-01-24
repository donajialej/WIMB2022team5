function Param_Est_SVIRD_uncert

clear all
close all
clc
format long 
warning('off','all')

global N alpha p gammasd gammasr gammavr gammavd m a b Y0

a = 1;
b = 90;
m = 10;

alpha = .75;
gammasd = .012*1/5;
gammavd = .0012*1/4;
gammasr = (1-.012)*1/5;
gammavr = (1-.0012)*1/4;
p = (.67-.56)/(b-a);

N = 10620000;

load Synth.txt

tdata = Synth(:,1);
Idata = Synth(:,2);
n = length(tdata);

disp('RECOVERED PARAMETERS FROM DETERMINISTIC METHOD - WITH RELATIVE DISCREPANCY - PREDICTOR-CORRECTOR')
disp('________________________________________________________________________________________________________________')
disp('j....... t1........ t2.........t3.........t4........ t5.........t6........t7.........t8.........t9.........t10..')
disp('________________________________________________________________________________________________________________')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_model(1)=0.26; theta_model(2)=-0.13; theta_model(3)=0.093; theta_model(4)=-0.036; theta_model(5)=0.05;
theta_model(6)=-0.12; theta_model(7)=0.14; theta_model(8)=-0.098; theta_model(9)=0.031; theta_model(10)=0.023;

vac_perc = .35;
inf_perc = .0003;
Y0(1) = (1 - vac_perc)*(1-inf_perc)*N;
Y0(2) = vac_perc*(1-inf_perc)*N;
I0 = inf_perc*N;

Y0(3) = I0*(1-vac_perc);
Y0(4) = I0*vac_perc;
Y0(5) = 0.012*Y0(3)+0.0012*Y0(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncertainty quantification - four parameter values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmb = 0;
NumCurves = 100;

curves = []; curves_results = [];  bt = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:NumCurves
theta = zeros(m,1);
theta(1) = .25;
theta(2) = -.1;
% theta(1:2) = theta_model(1:2);

nmb = nmb + 1;
    yirData = zeros(n,1);
%     yirData(1) = yi(1);
    
    for t = 1:n
        tau = abs(Idata(t));
        yirData(t,1) = poissrnd(tau,1,1);
    end
    
    curves = [curves (yirData)];

    ExpIncData = curves(1:n,iter);
    
    
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','InitDamping',1e5,'Disp','iter','TolX',10^(-2),'TolFun',10^(-2));
theta = lsqcurvefit(@(theta,tdata) Idata_function(theta,tdata),theta,tdata,ExpIncData,[],[],options);                               
                         
for j = 1:m
    Phatss(iter,j) = theta(j);
end    

%CREATION OF MATRIX FOR OUTPUT OF TABLE VALUES 
fprintf('%3d...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f\n',...
    nmb, theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7), theta(8), theta(9), theta(10)); 

curves_results = [curves_results (Idata_function(theta,tdata))];
bt = [bt (beta(tdata,theta))];

end

figure10 = figure; 
axes2 = axes('Parent',figure10,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all')

param_1=[mean(Phatss(:,1)) plims(Phatss(:,1),0.025) plims(Phatss(:,1),0.975)];

param_2=[mean(Phatss(:,2)) plims(Phatss(:,2),0.025) plims(Phatss(:,2),0.975)];

param_3=[mean(Phatss(:,3)) plims(Phatss(:,3),0.025) plims(Phatss(:,3),0.975)];

param_4=[mean(Phatss(:,4)) plims(Phatss(:,4),0.025) plims(Phatss(:,4),0.975)];

param_5=[mean(Phatss(:,5)) plims(Phatss(:,5),0.025) plims(Phatss(:,5),0.975)];

param_6=[mean(Phatss(:,6)) plims(Phatss(:,6),0.025) plims(Phatss(:,6),0.975)];

param_7=[mean(Phatss(:,7)) plims(Phatss(:,7),0.025) plims(Phatss(:,7),0.975)];

param_8=[mean(Phatss(:,8)) plims(Phatss(:,8),0.025) plims(Phatss(:,8),0.975)];

param_9=[mean(Phatss(:,9)) plims(Phatss(:,9),0.025) plims(Phatss(:,9),0.975)];

param_10=[mean(Phatss(:,10)) plims(Phatss(:,10),0.025) plims(Phatss(:,10),0.975)];


cad1=strcat('\theta_1=',num2str(param_1(end,1),2),'(95%CI:[',num2str(param_1(end,2),2),',',num2str(param_1(end,3),2),'])');
cad2=strcat('\theta_2=',num2str(param_2(end,1),2),'(95%CI:[',num2str(param_2(end,2),2),',',num2str(param_2(end,3),2),'])');
cad3=strcat('\theta_3=',num2str(param_3(end,1),2),'(95%CI:[',num2str(param_3(end,2),2),',',num2str(param_3(end,3),2),'])');
cad4=strcat('\theta_4=',num2str(param_4(end,1),2),'(95%CI:[',num2str(param_4(end,2),2),',',num2str(param_4(end,3),2),'])');
cad5=strcat('\theta_5=',num2str(param_5(end,1),2),'(95%CI:[',num2str(param_5(end,2),2),',',num2str(param_5(end,3),2),'])');

cad6=strcat('\theta_6=',num2str(param_6(end,1),2),'(95%CI:[',num2str(param_6(end,2),2),',',num2str(param_6(end,3),2),'])');
cad7=strcat('\theta_7=',num2str(param_7(end,1),2),'(95%CI:[',num2str(param_7(end,2),2),',',num2str(param_7(end,3),2),'])');
cad8=strcat('\theta_8=',num2str(param_8(end,1),2),'(95%CI:[',num2str(param_8(end,2),2),',',num2str(param_8(end,3),2),'])');
cad9=strcat('\theta_9=',num2str(param_9(end,1),2),'(95%CI:[',num2str(param_9(end,2),2),',',num2str(param_9(end,3),2),'])');
cad10=strcat('\theta_{10}=',num2str(param_10(end,1),2),'(95%CI:[',num2str(param_10(end,2),2),',',num2str(param_10(end,3),2),'])');

% suptitle(strcat('\fontsize{12}',cad1,';',cad2,';',cad3,';',cad4, ';',cad5));

subplot(2,5,1)
histogram(Phatss(:,1),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_1')
subplot(2,5,2)
histogram(Phatss(:,2),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_2')
subplot(2,5,3)
histogram(Phatss(:,3),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_3')
subplot(2,5,4)
histogram(Phatss(:,4),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_4')
subplot(2,5,5)
histogram(Phatss(:,5),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_5')

% suptitle(strcat('\fontsize{12}',cad6,';',cad7,';',cad8,';',cad9, ';',cad10));
subplot(2,5,6)
histogram(Phatss(:,6),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_6')
subplot(2,5,7)
histogram(Phatss(:,7),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_7')
subplot(2,5,8)
histogram(Phatss(:,8),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_8')
subplot(2,5,9)
histogram(Phatss(:,9),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_9')
subplot(2,5,10)
histogram(Phatss(:,10),'FaceColor',	[0.75, 0.75, 0])
xlabel('\fontsize{16}\theta_{10}')

% text1 = strcat('\fontsize{13}',cad1,{';  '},cad2,{';  '},cad3, {';  '},cad4,  {';  '},cad5);
% text(-4.0, 3.4, text1);
% 
% text2 = strcat('\fontsize{13}',cad6,{';  '},cad7,{';  '},cad8, {';  '},cad9,  {';  '},cad10);
% text(-4.0, -0.7, text2);


figure(figure10) 

disp([cad1 cad2 cad3 cad4 cad5])
disp([cad6 cad7 cad8 cad9 cad10])
tcalendar = datetime(2021,11,17) + caldays(1:length(tdata));

figure1 = figure;

line1 = plot(tcalendar,curves_results,'Color', 'b');

hold on 

% line2 = plot(tdata,curves,'Color',[0.9290, 0.6940, 0.1250]);

line3 = plot(tcalendar,Idata,'ro');
set(line3,'LineWidth',2)
hold on

line4 = plot(tcalendar,mean(curves_results,2),'-k');
set(line4,'LineWidth',2)

% line5 = line([20 20], [0 10000],'Color',[0.6350, 0.0780, 0.1840]);
% set(line5,'LineWidth',2)
% 
% text1 = strcat('\fontsize{13}',cad1,{';  '},cad2,{';  '},cad3, {';  '},cad4);
% text(9, 1750, text1);
% 
% axis([0 100 0 1600])

% ylim([0 1600])
xlim([datetime(2021,11,17) datetime(2022,02,17)])

legend([line1(1) line3(1) line4(1)],{ 'Reconstructed Incidence Curves',...
     'Reported Incidence Data', 'Mean of Reconstructed Incidence Curves',},'FontSize',14,'Location','best');

% title(strcat('\fontsize{14}','\it Forecasting from 40 Weeks of Incidence Data'))

% title(strcat('\fontsize{14}',cad1,';  ',cad2,';  ',cad3,';  ',cad4))


% title(strcat('\fontsize{12}',cad1,';',cad2,';',cad3,';',cad4, ';',cad5))
% title(strcat('\fontsize{12}',cad4, ';   ',cad5))

xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}Case incidence')

figure(figure1) 

figure6 = figure;
% line10 = fplot(@(t)beta(t,m,theta,a,b),[tdata(1) b],'-.^b');
line333 = plot(tcalendar,beta(tdata,theta_model),'g');
set(line333,'Linewidth',2)
hold on

line10 = plot(tcalendar, bt,'m');
hold on
line2a = plot(tcalendar,mean(bt,2),'-k');
set(line2a,'LineWidth',2)
% axis([0 90 0.0 1.0])

ylim([0.0 1.0])
xlim([datetime(2021,11,17) datetime(2022,02,17)])

legend([line333(1) line10(1) line2a(1)],{'Model Transmission Rate', 'Reconstructed Transmission Rates',...
    'Mean of Reconstructed Transmission Rates'},'FontSize',14,'Location','best');
% legend('Effective Reproduction Number','Initial Guess','best');
xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}\beta(t)')
figure(figure6)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FUNCTION DEFINITIONS

function dydt = svird(t,y,theta)
global N alpha p gammasd gammasr gammavr gammavd 
dydt = zeros(5,1);
dydt(1) = beta(t,theta).*y(1).*(y(3)+y(4))./(N - y(5)) - p.*y(1);
dydt(2) = p.*y(1) - (1-alpha)*beta(t,theta).*y(2).*(y(3)+y(4))./(N - y(5));
dydt(3) = beta(t,theta).*y(1).*(y(3)+y(4))./(N - y(5)) - (gammasr + gammasd)*y(3);
dydt(4) = (1-alpha).*beta(t,theta).*y(2).*(y(3)+y(4))./(N - y(5)) - (gammavr + gammavd)*y(4);
dydt(5) = gammasd.*y(3) + gammavd.*y(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yy = beta(t,theta) 
global m 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = Idata_function(theta,tdata)
global Y0 alpha N

[~,Y] = ode23s(@(t,y) svird(t,y,theta),tdata,Y0);

size(tdata)
size(beta(tdata,theta))
size(Y(:,1))

q = beta(tdata,theta).*Y(:,1).*(Y(:,3)+Y(:,4))./(N - Y(:,5)) +...
    (1-alpha)*beta(tdata,theta).*Y(:,2).*(Y(:,3)+Y(:,4))./(N - Y(:,5));

size(q)

end



function Param_Est_SVIRD_PCM_GA_normalized3
% A modified version with B = B(theta,u), GA data 7/9/21-11/25/21

clear all
close all
clc
format long 
warning('off','all')


global N alpha p gammasd gammasr gammavr gammavd delta1 delta2 m a b Y0 n K tdata scale1 scale2 scale3 rep_rate_D rep_rate_C

a = 1;
b = 140;

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

load GA_inc_cases_deaths_7_9_21_11_25_21_CDC_7day_average.txt
savename_results = 'fits/phat_beta_data_estimates_GA.mat';
tdata = GA_inc_cases_deaths_7_9_21_11_25_21_CDC_7day_average(:,1);
Idata = GA_inc_cases_deaths_7_9_21_11_25_21_CDC_7day_average(:,2);
Ddata = GA_inc_cases_deaths_7_9_21_11_25_21_CDC_7day_average(:,3);

K = length(tdata);

disp('RECOVERED EPIDEMIOLOGICAL PARAMETERS FROM PREDICTOR-CORRECTOR')
disp('________________________________________________________________________________________________________________')
disp('j....... t1........ t2.........t3.........t4........ t5.........t6........t7.........t8.........t9.........t10..')
disp('________________________________________________________________________________________________________________')

% Y0(1) = (N - 3942002 - 3580 - (1139396-21486))/N;
Y0(1) = (N - 3942002 - 3580)/N;
Y0(2) = 3942002/N;

Y0(3) = (3580 - 731)/N;
Y0(4) = 731/N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncertainty quantification 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%
tic;
%nmb = 0;
NumCurves = 100;
num_it = 130;
% 130
lambda0 = 1e0; 
% 0
power = 10; 

m = 10; % number of base functions for beta (epidemiological parameter)
n = 6;  % number of base functions for u (state variable) 
% 6
scale1 = 1e3;
% 1e3
scale2 = .05;
% 0.05
scale3 = 1;
% 1
rep_rate_D = .9;
% .9
rep_rate_C = .2;
% .2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predifine the size of vectors Alexandra was using
curvesI = zeros(K,NumCurves); curvesD = zeros(K,NumCurves); curves_results = zeros(K,NumCurves);
bt = zeros(K,NumCurves); rp = zeros(K,NumCurves);
SC = zeros(K,NumCurves); VC = zeros(K,NumCurves); ISC = zeros(K,NumCurves);  IVC = zeros(K,NumCurves); DC = zeros(K,NumCurves);
Phatss = zeros(NumCurves,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate the random numbers first for reproducibility:
rng(0) % fixing the seed for reproducibility
tau1 = repmat(abs(Idata),[1,NumCurves]);
yirIdata = poissrnd(tau1);
tau2 = repmat(abs(Ddata),[1,NumCurves]);
yirDdata = poissrnd(tau2);


parfor iter = 1:NumCurves
    %nmb = nmb + 1;
    %yirIdata = zeros(K,1);
    %yirDdata = zeros(K,1);
    %tau1 = abs(Idata);
    %yirIdata = poissrnd(tau1);
    %tau2 = abs(Ddata);
    %yirDdata = poissrnd(tau2);

    curvesI(:,iter) = yirIdata(:,iter);
    curvesD(:,iter) = yirDdata(:,iter);

    ExpIncData = curvesI(:,iter);
    ExpCumDeaths = cumsum(curvesD(:,iter));
    ExpIncDeaths = curvesD(:,iter);
    
    % Normalize the data
    Data = zeros(2*K,1);
    DataN = zeros(2*K,1);
    Data(1:K,1) = ExpIncData;
    DataN(1:K,1) = ExpIncData./N;
    Data(K+1:2*K,1) = scale1*ExpIncDeaths;
    DataN(K+1:2*K,1) = scale1*ExpIncDeaths./N;
    ExpCumDeathsN = ExpCumDeaths/N;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    theta0 = zeros(m,1);
    theta0(1,1) = .5; 
    theta1 = theta0;
    u = u0f(theta0,ExpCumDeathsN,n,a,b,Y0,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options1 = optimoptions('lsqnonlin','MaxIterations',0,...
        'Algorithm','levenberg-marquardt','InitDamping',1e5,'Disp','off');
    
    options2 = optimoptions('lsqcurvefit','MaxIterations',0,...
        'Algorithm','levenberg-marquardt','InitDamping',1e5,'Disp','off'); 
    
    for lp = 1:num_it  
        % Evalaute all four Jacobians
        [~,~,~,~,~,~,jacobian] = lsqnonlin(@(theta) G(theta,u,ExpCumDeathsN,m,a,b,n,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, K, tdata, scale2),theta1,[],[],options1);

        G_prime_theta = jacobian;
    
        [~,~,~,~,~,~,jacobian] = lsqcurvefit(@(theta,tdata) B(theta,u,ExpCumDeathsN,m,a,b,n,tdata, alpha, gammasd, gammavd, K, scale1, rep_rate_D, rep_rate_C),theta1,tdata,DataN,[],[],options2);
     
        B_prime_theta = jacobian;

        lambda = lambda0/lp^power;
        % lambda = lambda0/exp(factor*lp);
    
        step_theta = -(G_prime_theta'*G_prime_theta + B_prime_theta'*B_prime_theta + lambda*eye(m))\...
            (G_prime_theta'*G(theta1,u,ExpCumDeathsN,m,a,b,n,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, K, tdata, scale2)+ scale3*B_prime_theta'*(B(theta1,u,ExpCumDeathsN,m,a,b,n,tdata, alpha, gammasd, gammavd, K, scale1, rep_rate_D, rep_rate_C) - DataN) + lambda*(theta1 - theta0));
        theta1 = step_theta + theta1;

        [~,~,~,~,~,~,jacobian] = lsqnonlin(@(u) G(theta1,u,ExpCumDeathsN,m,a,b,n,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, K, tdata, scale2),u,[],[],options1);

        G_prime_u = jacobian;
    
        [~,~,~,~,~,~,jacobian] = lsqcurvefit(@(u,tdata) B(theta1,u,ExpCumDeathsN,m,a,b,n,tdata, alpha, gammasd, gammavd, K, scale1, rep_rate_D, rep_rate_C),u,tdata,DataN,[],[],options2);
    
        B_prime_u = jacobian;
    
        step_u = -(G_prime_u'*G_prime_u + B_prime_u'*B_prime_u)\(G_prime_u'*G(theta1,u,ExpCumDeathsN,m,a,b,n,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, K, tdata, scale2)...
            + scale3*B_prime_u'*(B(theta1,u,ExpCumDeathsN,m,a,b,n,tdata, alpha, gammasd, gammavd, K, scale1, rep_rate_D, rep_rate_C) - DataN));
        u = .1*step_u + u;

    end
 
    
    Phatss(iter,:) = theta1;
       

    %CREATION OF MATRIX FOR OUTPUT OF TABLE VALUES 
    fprintf('%3d...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f\n',...
        iter, theta1(1), theta1(2), theta1(3), theta1(4), theta1(5), theta1(6), theta1(7), theta1(8), theta1(9), theta1(10)); 

    yd = round(rep_rate_D*N.*(gammasd.*y(3,tdata,u,n,a,b) + gammavd.*y(4,tdata,u,n,a,b)));

    yc = round(rep_rate_C*N.*(beta(tdata,theta1,m,a,b).*y(1,tdata,u,n,a,b).*(y(3,tdata,u,n,a,b) + y(4,tdata,u,n,a,b))./(1-ExpCumDeathsN)...
          + (1 - alpha).*beta(tdata,theta1,m,a,b).*y(2,tdata,u,n,a,b).*(y(3,tdata,u,n,a,b) + y(4,tdata,u,n,a,b))./(1-ExpCumDeathsN)));
      
    repr = beta(tdata,theta1,m,a,b).*y(1,tdata,u,n,a,b)./((1-ExpCumDeathsN).*(gammasr + gammasd))...
          + (1 - alpha).*beta(tdata,theta1,m,a,b).*y(2,tdata,u,n,a,b)./((1-ExpCumDeathsN).*(gammavr + gammavd));     

    y1 = round(N.*y(1,tdata,u,n,a,b));
    y2 = round(N.*y(2,tdata,u,n,a,b));
    y3 = round(N.*y(3,tdata,u,n,a,b));
    y4 = round(N.*y(4,tdata,u,n,a,b));

    curves_results(:,iter) = yc;
    bt(:,iter) = beta(tdata,theta1,m,a,b);
    rp(:,iter) = repr;

    SC(:,iter) = y1;
    VC(:,iter) = y2;
    ISC(:,iter) = y3;
    IVC(:,iter) = y4;
    DC(:,iter) = yd;

end

tcalendar = datetime(2021,07,9) + caldays(1:length(tdata));

% Save simulation
save(savename_results,'Phatss','tcalendar','bt','rp','SC','VC','ISC','IVC','DC','curves_results','Idata','Ddata')


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

disp([cad1 cad2 cad3 cad4 cad5])
disp([cad6 cad7 cad8 cad9 cad10])

% Figures

%%%%%% one figure
figure10 = figure;%('PaperPosition',[.25 .25 11.7 4]);
axes2 = axes('Parent',figure10,... 
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
box(axes2,'on');
hold(axes2,'all')
%set(gca,'fontsize', 16);

subplot(4,5,1)
histogram(Phatss(:,1),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_1')
subplot(4,5,2)
histogram(Phatss(:,2),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_2')
subplot(4,5,3)
histogram(Phatss(:,3),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_3')
subplot(4,5,4)
histogram(Phatss(:,4),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_4')
subplot(4,5,5)
histogram(Phatss(:,5),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_5')

subplot(4,5,6)
histogram(Phatss(:,6),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_6')
subplot(4,5,7)
histogram(Phatss(:,7),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_7')
subplot(4,5,8)
histogram(Phatss(:,8),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_8')
subplot(4,5,9)
histogram(Phatss(:,9),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_9')
subplot(4,5,10)
histogram(Phatss(:,10),'FaceColor',	[0.75, 0.75, 0])
xlabel('\theta_{10}')


subplot(4,5,[11 12 13 14 15 16 17 18 19 20])
line10 = plot(tcalendar, bt,'m');
hold on
line2a = plot(tcalendar,mean(bt,2),'-^k');
line5x = line([datetime(2021,07,9) datetime(2021,11,25)], [1 1],'Color',[0.6350, 0.0780, 0.1840]);
line10x = plot(tcalendar, rp,'g');
line2x = plot(tcalendar,mean(rp,2),'-^k');

set(line10,'LineWidth',1)
set(line2a,'LineWidth',1)
set(line5x,'LineWidth',2)
set(line10x,'LineWidth',1)
set(line2x,'LineWidth',1)


ylim([0.0 2.0])
xlim([datetime(2021,07,9) datetime(2021,11,25)])

legend([line10(1) line2a(1) line10x(1) line2x(1)],{'Reconstructed \beta(t)','Mean Reconstructed \beta(t)',... 
    'Effective Reproduction Number', 'Mean of Effective Reproduction Number'},'FontSize',18,'Location','best');

xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}\beta(t)')

%set(gca,'fontsize', 16);
figure(figure10) 

saveas(gcf,'figures/Figure8','epsc')
savefig('figures/Figure8')


toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FUNCTION DEFINITIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = G(theta,u,D,m,a,b,n,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, K, tdata, scale2)
%global alpha p gammasd gammasr gammavr gammavd delta1 delta2 K tdata scale2

yR = 1 - D;

for ii = 1:4
yR = yR - y(ii,tdata,u,n,a,b);
end
for jj = 1:K
if yR(jj) < 0; yR(jj) = 0; end
end

g(1:K,1) = dydt(1,tdata,u,n,a,b) - (-beta(tdata,theta,m,a,b).*y(1,tdata,u,n,a,b).*(y(3,tdata,u,n,a,b)+y(4,tdata,u,n,a,b))./(1-D)...
    - p.*y(1,tdata,u,n,a,b) + delta1*yR + delta2*y(2,tdata,u,n,a,b));

g(K+1:2*K,1) = dydt(2,tdata,u,n,a,b) - (p.*y(1,tdata,u,n,a,b) - (1-alpha)*beta(tdata,theta,m,a,b).*y(2,tdata,u,n,a,b).*(y(3,tdata,u,n,a,b)...
    + y(4,tdata,u,n,a,b))./(1-D) - delta2*y(2,tdata,u,n,a,b)); 

g(2*K+1:3*K,1) = dydt(3,tdata,u,n,a,b) - (beta(tdata,theta,m,a,b).*y(1,tdata,u,n,a,b).*(y(3,tdata,u,n,a,b)+y(4,tdata,u,n,a,b))./(1-D)...
     - (gammasr + gammasd)*y(3,tdata,u,n,a,b));

g(3*K+1:4*K,1) = dydt(4,tdata,u,n,a,b) - ((1-alpha).*beta(tdata,theta,m,a,b).*y(2,tdata,u,n,a,b).*(y(3,tdata,u,n,a,b)+y(4,tdata,u,n,a,b))./(1-D)...
     - (gammavr + gammavd)*y(4,tdata,u,n,a,b));

g = scale2*g;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bb = B(theta,u,D,m,a,b,n,tdata, alpha, gammasd, gammavd, K, scale1, rep_rate_D, rep_rate_C)
%global  tdata alpha gammasd gammavd K scale1 rep_rate_D rep_rate_C
    
  bb(1:K,1) =  rep_rate_C*((y(1,tdata,u,n,a,b) + y(2,tdata,u,n,a,b).*(1 - alpha)).*beta(tdata,theta,m,a,b).*(y(3,tdata,u,n,a,b) + y(4,tdata,u,n,a,b))./(1-D));
  bb(K+1:2*K,1) = rep_rate_D*scale1*(gammasd.*y(3,tdata,u,n,a,b) + gammavd.*y(4,tdata,u,n,a,b));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yy = y(i,t,u,n,a,b)
%global n a b 

yy = 0;
   for j = 1:n
       yy = yy + u(j+(i-1)*n,1).*leg(j-1,t,a,b);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yy = dydt(i,t,u,n,a,b)
%global n a b 

yy = 0;
   for j = 1:n
       yy = yy + u(j+(i-1)*n).*leg_prime(j-1,t,a,b);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yy = beta(t,theta,m,a,b) 
%global m a b 

       yy = 0;
   for j = 1:m
       yy = yy + theta(j).*leg(j-1,t,a,b);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = leg(j,t,A,B)

x = (2.*t - A - B)./(B - A);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function P_prime = leg_prime(j,t,A,B)

x = (2.*t - A - B)./(B - A);

if j == 0
    P1_prime = 0; P_prime = P1_prime;
elseif j == 1 
    P2_prime = 2/(B-A); P_prime = P2_prime;
else   
   P1 = 1; P2 = x; P2_prime = 2/(B-A);
       for k = 2:j
           P3 = ((2*(k-1)+1).*x.*P2 - (k-1).*P1)./k;
           P3_prime = 2*k*P2/(B-A) + x.*P2_prime;
           P1 = P2; P2 = P3; 
           P2_prime = P3_prime;
       end
           P_prime = P3_prime;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = svird(t,y,theta,D, alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m,a,b)
%global alpha p gammasd gammasr gammavr gammavd delta1 delta2 tdata

Dsp = spline(tdata,D,t);

dy = zeros(4,1);

y5 = 1 - sum(y(1:4)) - Dsp;

dy(1) = -beta(t,theta,m,a,b).*y(1).*(y(3)+y(4))./(1-Dsp) - p.*y(1) + delta1*y5 + delta2*y(2);
dy(2) = p.*y(1) - (1-alpha)*beta(t,theta,m,a,b).*y(2).*(y(3)+y(4))./(1-Dsp) - delta2*y(2);
dy(3) = beta(t,theta,m,a,b).*y(1).*(y(3)+y(4))./(1-Dsp) - (gammasr + gammasd)*y(3);
dy(4) = (1-alpha).*beta(t,theta,m,a,b).*y(2).*(y(3)+y(4))./(1-Dsp) - (gammavr + gammavd)*y(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yy = u0f(theta,D,n,a,b,Y0,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m)
%global n

yy = zeros(4*n,1);

for i = 1:4
    for  j = 1:n
           yy(j+(i-1)*n,1) = (2*(j-1)+1)*integral(@(x)fun(x,j,i,theta,D,a,b,Y0,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m),-1,1)/2;
    end

end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun(x,j,i,theta,D,a,b,Y0,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m)
%global a b Y0 

t = ((b-a).*x+a+b)./2;

[T,Y] = ode23s(@(t,y) svird(t,y,theta,D,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m,a,b),[a b],Y0);

y = spline(T,Y(:,i),t).*leg(j-1,x,-1,1);

end
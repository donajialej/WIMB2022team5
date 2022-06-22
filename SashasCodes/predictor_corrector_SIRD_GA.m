function predictor_corrector_SIRD_GA
%close all 
clc  
format long 
warning('off','all')

% load SLcumcases.txt;
% load GUcumcases.txt;
% load LIBcumcases.txt;
% load SFcumcases.txt;
% load Mcumcases.txt;
%load Georgiacumcases.txt;

global Cdata Cdata_inc tdata n gamma nu a b N m NP
gamma = 1/5;
% nu = 0.005;
nu = 0.012;

% tdata = SLcumcases(:,1); % alpha0 = 1e-2; wK = 1e0; w = 1e9; K0 = 20000
% Cdata = SLcumcases(:,2); 
% tdata = GUcumcases(:,1); % alpha0 = 1e-3; wK = 1e0; w = 1e9; K0 = 7000
% Cdata = GUcumcases(:,2);
% tdata = LIBcumcases(:,1);
% Cdata = LIBcumcases(:,2);
% tdata = SFcumcases(:,1);
% Cdata = SFcumcases(:,2);
% tdata = Mcumcases(:,1); %  alpha0 = 1e-3; wK = 1e0; w = 1e9; K0 = 5000
% Cdata = Mcumcases(:,2);
tdata = Georgiacumcases(:,1);
Cdata = Georgiacumcases(:,2);
Cdata_inc = [Cdata(1,1);diff(Cdata(:,1))];

n = length(tdata);
m = 10; % number of Legendre polynomials
num_it = 13;
N = 15;
NP = 10620000; % population of Georgia was 10.62 million in 2019
%--------------------------------------------------------------------------
% DETERMENISTIC RECOVERY - PREDICTOR-CORRECTOR ALGORITHM
%--------------------------------------------------------------------------
% TO DISPLAY ITERATIONS FOR DETERMENISTIC METHOD 

disp('RECOVERED PARAMETERS FROM DETERMINISTIC METHOD - WITH RELATIVE DISCREPANCY - PREDICTOR-CORRECTOR')
disp('____________________________________________________________________________________________________________________________')
disp('j....... t1........ t2.........t3.........t4........ t5.........t6........t7.........t8.........t9.........t10.........RD...')
disp('____________________________________________________________________________________________________________________________')
 

a = tdata(1)-1; b = tdata(n); 

u = u0;

alpha0 = 1e-7;
% alpha0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncertainty quantification - four parameter values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmb = 0;
NumCurves = 3;

curves = []; curves_results = []; initial = []; bt = []; rp = [];
curves_results_cum = [];

curves_S = [];
curves_I = [];
curves_R = [];
curves_D = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:NumCurves
theta  = zeros(m,1);
theta(1) = .1 + (1 -.1).*rand(1,1);
theta0 = zeros(m,1);
theta0(1) = .1 + (1 -.1).*rand(1,1);

  nmb = nmb + 1;
    yirData = zeros(n,1);
%     yirData(1) = yi(1);
    
    for t = 1:n
        tau = abs(Cdata_inc(t));
        yirData(t,1) = poissrnd(tau,1,1);
    end
    
    curves = [curves (yirData)];
    
    init = theta(1)*ones(n,1);
    initial = [initial (init)];

    ExpIncData = curves(1:n,iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_vect = zeros(n,1); B_vect = zeros(n,1); g = zeros(n,1);
G_prime_theta = zeros(n,m); G_prime_u = zeros(n,2*N); B_prime_u = zeros(n,2*N); 

for i = 1:n
    G_vect(i,1) = G(theta,u,i);
    for k = 1:m
        G_prime_theta(i,k) = Gthki(u,k,i);
    end
end

for j = 1:num_it

% alpha = alpha0/j^.25;
alpha = alpha0;
% alpha = alpha0/j;
% alpha = alpha0/exp(j);
    
    step_theta = -(G_prime_theta'*G_prime_theta + alpha*eye(m))\(G_prime_theta'*(G_vect - g)...
        + alpha*eye(m)*(theta - theta0));
    theta = step_theta + theta;
  
%     for l=1:m
%     if (theta(l)<0)
%         theta(l) = theta0(l);
%     end
%     end
    
for i = 1:n
    G_vect(i,1) = G(theta,u,i);
    B_vect(i,1) = B(u,i);
    for k = 1:2*N
        B_prime_u(i,k) = Buki(k,i);
        G_prime_u(i,k) = Guki(theta,u,k,i);
    end
end
    
    step_u = -(G_prime_u'*G_prime_u + B_prime_u'*B_prime_u)\(G_prime_u'*(G_vect - g)...
        + B_prime_u'*(B_vect - ExpIncData));
    u = step_u + u;
    
for i = 1:n
    G_vect(i,1) = G(theta,u,i);
    B_vect(i,1) = B(u,i);
    for k = 1:m
        G_prime_theta(i,k) = Gthki(u,k,i);
    end
end
    
RD1 = norm(B_vect - ExpIncData)/norm(Cdata_inc);
RD2 = norm(G_vect - g)/norm(Cdata_inc);
RD = RD1 + RD2;

% if (RD > 1.01*RD0)
% if alpha < 1.25*1e-5
%         fprintf('Stopping Time!\n');
%         break;
% end

 RD0 = RD;
 
end

% Phatss is the matrix with the estimated parameters
  
for j = 1:m
    Phatss(iter,j) = theta(j);
end    

%CREATION OF MATRIX FOR OUTPUT OF TABLE VALUES 
fprintf('%3d...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f\n',...
    nmb, theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7), theta(8), theta(9), theta(10),RD); 

% Matrix with the estimated values of the reproduction number

% for i=1:n
% y(iter,i) = reproduction(tdata(i),m,theta,a,b);
% end

curves_results = [curves_results (E(tdata,u))];
curves_results_cum = [curves_results_cum (C(tdata,u))];
bt = [bt (beta(tdata,m,theta,a,b))];
rp = [rp (reproduction(tdata,m,theta,u,a,b))];
curves_S = [curves_S (SF(tdata,u))];
curves_I = [curves_I (IF(tdata,u))];
curves_R = [curves_R (RF(tdata,u))];
curves_D = [curves_D (DF(tdata,u))];

end

RD1
RD2
alpha

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
tcalendar = datetime(2020,03,04) + caldays(1:length(tdata));

figure1 = figure;

line1 = plot(tcalendar,curves_results,'Color', 'g');

hold on 

% line2 = plot(tdata,curves,'Color',[0.9290, 0.6940, 0.1250]);

line3 = plot(tcalendar,Cdata_inc,'ro');
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

ylim([0 1600])
xlim([datetime(2020,03,04) datetime(2020,06,02)])

legend([line1(1) line3(1) line4(1)],{ 'Reconstructed Incidence Curves',...
     'Reported Incidence Data', 'Mean of Reconstructed Incidence Curves',},'FontSize',14,'Location','best');

% title(strcat('\fontsize{14}','\it Forecasting from 40 Weeks of Incidence Data'))

% title(strcat('\fontsize{14}',cad1,';  ',cad2,';  ',cad3,';  ',cad4))


% title(strcat('\fontsize{12}',cad1,';',cad2,';',cad3,';',cad4, ';',cad5))
% title(strcat('\fontsize{12}',cad4, ';   ',cad5))

xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}Case incidence')

figure(figure1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure2 = figure;

line1a = plot(tcalendar,curves_results_cum,'Color', 'c');

hold on 

% line2 = plot(tdata,curves,'Color',[0.9290, 0.6940, 0.1250]);

line3a = plot(tcalendar,Cdata,'ro');
set(line3,'LineWidth',2)
hold on

line4a = plot(tcalendar,mean(curves_results_cum,2),'-k');
set(line4a,'LineWidth',2)

% line5 = line([20 20], [0 10000],'Color',[0.6350, 0.0780, 0.1840]);
% set(line5,'LineWidth',2)
% 
% text1 = strcat('\fontsize{13}',cad1,{';  '},cad2,{';  '},cad3, {';  '},cad4);
% text(9, 1750, text1);
% 
% axis([0 100 0 1600])

% ylim([0 1600])
xlim([datetime(2020,03,04) datetime(2020,06,02)])

legend([line1a(1) line3a(1) line4a(1)],{ 'Reconstructed Cumulative Curves',...
     'Reported Cumulative Data', 'Mean of Reconstructed Cumulative Curves',},'FontSize',14,'Location','best');

% title(strcat('\fontsize{14}','\it Forecasting from 40 Weeks of Incidence Data'))

% title(strcat('\fontsize{14}',cad1,';  ',cad2,';  ',cad3,';  ',cad4))


% title(strcat('\fontsize{12}',cad1,';',cad2,';',cad3,';',cad4, ';',cad5))
% title(strcat('\fontsize{12}',cad4, ';   ',cad5))

xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}Cumulative Cases')

figure(figure2) 

figure6 = figure;
% line10 = fplot(@(t)beta(t,m,theta,a,b),[tdata(1) b],'-.^b');
line2 = plot(tcalendar,initial, 'Color', [0.8500, 0.3250, 0.0980]);
hold on
% fplot(@(t)r0(t),[a b],'-.dr');
% line10 = plot(tdata, bt,'-.^b');
line10 = plot(tcalendar, bt,'b');
hold on
line2a = plot(tcalendar,mean(bt,2),'-k');
set(line2a,'LineWidth',2)
% axis([0 90 0.0 1.0])

ylim([0.0 1.0])
xlim([datetime(2020,03,04) datetime(2020,06,02)])

legend([line2(1) line10(1) line2a(1)],{'Initial Approximations', 'Reconstructed Transmission Rates',...
    'Mean of Reconstructed Transmission Rates'},'FontSize',14,'Location','best');
% legend('Effective Reproduction Number','Initial Guess','best');
xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}\beta(t)')
figure(figure6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure66 = figure;

line100 = plot(tcalendar, rp,'m');
hold on
line22a = plot(tcalendar,mean(rp,2),'-.^k');
set(line22a,'LineWidth',2)
hold on
line222 = line([tcalendar(1) tcalendar(end)], [1 1], 'Color',[0.4660, 0.6740, 0.1880]);
set(line222,'LineWidth',2)

ylim([0.0 6.0])
xlim([datetime(2020,03,04) datetime(2020,06,02)])

legend([line100(1) line22a(1)],{'Reconstructed Effective Reproduction Numbers',...
    'Mean of Effective Reproduction Numbers'},'FontSize',18,'Location','best');
% legend('Effective Reproduction Number','Initial Guess','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14}\it{R}(t)')
figure(figure66)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure7 = figure;
subplot(2,2,1)
line77 = plot(tcalendar,curves_S, '-g');
hold on 
line77a = plot(tcalendar,mean(curves_S,2),'-k');
set(line77a,'LineWidth',2);
legend([line77(1) line77a(1)],{'Reconstructed S(t)', ...
    'Mean of S(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} S(t)')
% hold on 
subplot(2,2,2)
line78 = plot(tcalendar,curves_I, 'Color', [0.6350, 0.0780, 0.1840]);
hold on 
line78a = plot(tcalendar,mean(curves_I,2),'-k');
set(line78a,'LineWidth',2);
legend([line78(1) line78a(1)],{'Reconstructed I(t)', ...
    'Mean of  I(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} I(t)')
% hold on 
subplot(2,2,3)
line79 = plot(tcalendar,curves_R, '-m');
hold on 
line79a = plot(tcalendar,mean(curves_R,2),'-k');
set(line79a,'LineWidth',2);
legend([line79(1) line79a(1)],{'Reconstructed R(t)', ...
    'Mean of  R(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} R(t)')
% hold on 
subplot(2,2,4)
line76 = plot(tcalendar,curves_D, '-c');
hold on 
line76a = plot(tcalendar,mean(curves_D,2),'-k');
set(line76a,'LineWidth',2);
legend([line76(1) line76a(1)],{'Reconstructed D(t)', ...
    'Mean of D(t)'},'FontSize',14,'Location','best');
xlabel('\fontsize{14}Time (days)');
ylabel('\fontsize{14} D(t)')
figure(figure7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       function y = dIdt(t,u)
       global a b N
       
       Af = u(1:N); Bf = u(N+1:2*N);
       
       y = 0;
       for j = 1:N
           
         y = y + Af(j)*(cos(2*pi*j*(t-a)./(b-a))-1) + Bf(j)*sin(2*pi*j*(t-a)./(b-a));
         
       end
       
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       function y = I(t,u)
       global a b N
       
       Af = u(1:N); Bf = u(N+1:2*N);
       
       y = 0;
       for j = 1:N
         y = y + Af(j)*((b-a)*sin(2*pi*j*(t-a)./(b-a))/(2*pi*j) - (t-a))...
               - Bf(j)*(b-a)/(2*pi*j)*(cos(2*pi*j*(t-a)./(b-a))-1);
       end
       
       end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       function y = E(t,u)
       global gamma
       
       y = dIdt(t,u) + gamma*I(t,u);
       
       end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       function y = IntI(t,u)
       global a b N
       
       Af = u(1:N); Bf = u(N+1:2*N);
       
       y = 0;
       for j = 1:N
           y = y + Af(j)*(-(b-a)^2*(cos(2*pi*j*(t-a)./(b-a)) - 1)./(2*pi*j)^2 - (t-a).^2./2) ...
               - Bf(j)*(b-a)/(2*pi*j)*((b-a)/(2*pi*j)*sin(2*pi*j*(t-a)./(b-a))-(t-a));
           
           
       end
       
       end
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = C(t,u)
       global gamma
       
       y = I(t,u) + gamma*IntI(t,u);
       
       end  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = SF(t,u)
       global NP gamma
       
       y = NP - I(t,u)- gamma*IntI(t,u);
%        y = y/NP;
       
       end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = IF(t,u)
       global NP 
       
%        y = I(t,u)/NP;
         y = I(t,u);
       
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = RF(t,u)
       global NP gamma nu
       
       y = (1-nu)*gamma*IntI(t,u);
%        y = y/NP;
       
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = DF(t,u)
       global NP gamma nu
       
       y = nu*gamma*IntI(t,u);
%        y = y/NP;
       
       end                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coef = u0
global a b N

coef = zeros(2*N,1);
Af = zeros(N,1); Bf = zeros(N,1);

for j = 1:N
    
Af(j,1) = 2*integral(@(x)fun1(x,j),a,b)/(b-a);
Bf(j,1) = 2*integral(@(x)fun2(x,j),a,b)/(b-a);

end

coef(1:N) = Af; coef(N+1:2*N) = Bf;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun1(x,j)
global a b
       y = initial_guess_dIdt(x).*cos(2*pi*j.*(x-a)./(b-a));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun2(x,j)
global a b
       y = initial_guess_dIdt(x).*sin(2*pi*j.*(x-a)./(b-a));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = initial_guess_I(t)
global gamma Cdata_inc tdata n

yy = zeros(n,1);
for k = 1:n
    yy(k) = exp(-gamma*tdata(k)).*sum(exp(gamma*tdata(1:k)).*Cdata_inc(1:k));
end

pp = spline(tdata,yy);
y = ppval(pp,t);

end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = initial_guess_dIdt(t)
global gamma Cdata_inc tdata

pp = spline(tdata,Cdata_inc);
yy = ppval(pp,t);

        y = yy - gamma*initial_guess_I(t);

end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = initial_guess_IntI(t)
global Cdata gamma tdata

%        y = integral(@(x)initial_guess_I(x),a,t);
         pp = spline(tdata,Cdata);
         yy = ppval(pp,t);
         y = (yy - initial_guess_I(t))./gamma;

end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = initial_guess_cumulative(t)
global gamma

    y = initial_guess_I(t) + gamma*initial_guess_IntI(t);


end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = initial_guess_incidence(t)
global gamma

    y = initial_guess_dIdt(t) + gamma*initial_guess_I(t);

end 

% FUNCTIONS FOR THE ITERATIVE PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = B(u,i)

global gamma tdata

y = dIdt(tdata(i),u) + gamma*I(tdata(i),u);
 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Buki(k,i)

global gamma tdata

y = dIdt_u(tdata(i),k) + gamma*I_u(tdata(i),k);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = G(theta,u,i)

global NP a b m tdata gamma nu 

y = B(u,i) - beta(tdata(i),m,theta,a,b).*I(tdata(i),u).*(NP...
    - I(tdata(i),u) - gamma*IntI(tdata(i),u))./(NP - nu*gamma*IntI(tdata(i),u));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = Guki(theta,u,k,i) 

global tdata NP a b nu gamma m

y = Buki(k,i) - beta(tdata(i),m,theta,a,b).*(I(tdata(i),u).*((- I_u(tdata(i),k) - gamma*IntI_u(tdata(i),k)).*(NP...
- nu*gamma*IntI(tdata(i),u))-(NP - I(tdata(i),u) - gamma*IntI(tdata(i),u)).*(- nu*gamma*IntI_u(tdata(i),k)))...
    ./(NP - nu*gamma*dIdt(tdata(i),u)).^2 + I_u(tdata(i),k).*(NP...
    - I(tdata(i),u) - gamma*IntI(tdata(i),u))./(NP - nu*gamma*IntI(tdata(i),u)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = dIdt_u(t,k)

global a b N

if (1<=k)&&(k<=N) 
    y = cos(2*pi*k*(t-a)./(b-a))-1; 
elseif (N+1<=k)&&(k<=2*N)
    y = sin(2*pi*(k-N)*(t-a)./(b-a));    
end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = I_u(t,k)   
global a b N

if (1<=k)&&(k<=N) 
    y = (b-a)*sin(2*pi*k*(t-a)./(b-a))/(2*pi*k) - (t-a); 
elseif (N+1<=k)&&(k<=2*N)
    y = -(b-a)/(2*pi*k)*(cos(2*pi*k*(t-a)./(b-a))-1);   
end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = IntI_u(t,k)  

global a b N

if (1<=k)&&(k<=N) 
    y =  -(b-a)^2*(cos(2*pi*k*(t-a)./(b-a)) - 1)./(2*pi*k)^2 - (t-a).^2./2;
elseif (N+1<=k)&&(k<=2*N)
    y = -(b-a)/(2*pi*k)*((b-a)/(2*pi*k)*sin(2*pi*k*(t-a)./(b-a))-(t-a));   
end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Gthki(u,k,i)
global NP tdata gamma nu
      
  bth = bthki(i);
  y = - bth(k).*I(tdata(i),u).*(NP - I(tdata(i),u) - gamma*IntI(tdata(i),u))./(NP - nu*gamma*IntI(tdata(i),u));
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = bthki(i)
global tdata m a b
y = zeros(m,1);     
       
   for j = 1:m
       y(j) = leg(j-1,tdata(i),a,b);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = beta(t,m,theta,a,b)     
       y = 0;
   for j = 1:m
       y = y + theta(j).*leg(j-1,t,a,b);
   end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = reproduction(t,m,theta,u,a,b)    

global gamma NP
   
y = beta(t,m,theta,a,b).*SF(t,u)./((NP - DF(t,u))*gamma);
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = horiz(t)  

y = 0*t + 1;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = leg(j,t,a,b)
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

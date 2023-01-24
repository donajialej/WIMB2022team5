function Param_Est_SVIRD_PCM_CAL_mult_rep_rates
% A modified version with B = B(theta,u), GA data 7/9/21-11/25/21
tic;
clear all
close all
clc
format long 
warning('off','all')


global N alpha p gammasd gammasr gammavr gammavd delta1 delta2 m a b Y0 K tdata scale1 scale2 scale3 rep_rate_D 

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

N = 39237836; % CA population

load California_7_9_21_11_25_21.txt

tdata = California_7_9_21_11_25_21(:,1);
Idata = California_7_9_21_11_25_21(:,2);
Ddata = California_7_9_21_11_25_21(:,3);
Dcum = cumsum(Ddata(1:length(tdata)));

K = length(tdata);

disp('RECOVERED EPIDEMIOLOGICAL PARAMETERS FROM PREDICTOR-CORRECTOR')
disp('________________________________________________________________________________________________________________')
disp('j....... t1........ t2.........t3.........t4........ t5.........t6........t7.........t8.........t9.........t10..')
disp('________________________________________________________________________________________________________________')

Y0(1) = (N - 20086693 - 25039)/N;
Y0(2) = 20086693/N;

Y0(3) = (25039 - 5116)/N;
Y0(4) = 5116/N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%

%nmb = 0;
NumCurves = 10;
% num_it = 58;
% 41
% lambda0 = 1e3; 
% 1e2
power = 7; 

m = 10; % number of base functions for beta (epidemiological parameter)
% n = 12;  % number of base functions for u (state variable) 
% 10
scale1 = 1e3;
% 1e3
scale2 = .1;
% 0.05
scale3 = 1;
% 1
rep_rate_D = .9;
% .9
rep_rate = .1;
% .2

% Normalize the data
    DataN = zeros(2*K,1);
    DataN(1:K,1) = Idata./N;
    DataN(K+1:2*K,1) = scale1*Ddata./N;
    ExpCumDeathsN = Dcum/N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rp = zeros(K,NumCurves);

for iter = 1:NumCurves
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      if iter==1
           num_it=41; lambda0=1e3; rep_rate_C=.1;n=12;
       elseif iter==2
           num_it=58; lambda0=1e3;rep_rate_C=.2;n=12;
       elseif iter==3
           num_it=45; lambda0=1e2;rep_rate_C=.3;n=11;
       elseif iter==4
           num_it=45; lambda0=1e2;rep_rate_C=.4;n=10;
       elseif iter==5
           num_it=40; lambda0=1e2;rep_rate_C=.5;n=10;
       elseif iter==6
           num_it=40; lambda0=1e2;rep_rate_C=.6;n=10;
       elseif iter==7
           num_it=45; lambda0=1e2;rep_rate_C=.7;n=10;
       elseif iter==8
           num_it=45; lambda0=1e2;rep_rate_C=.8;n=10;
       elseif iter==9 
           num_it=50; lambda0=1e2;rep_rate_C=.9;n=10;
       else 
        num_it=50; lambda0=1e2;rep_rate_C=.9;n=10;
       end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    theta0 = zeros(m,1);
    theta0(1,1) = .5; 
    theta1 = theta0;
    u = u0f(theta0,ExpCumDeathsN,n,a,b,Y0,alpha, p, gammasd, gammasr, gammavr, gammavd, delta1, delta2, tdata,m);
    
%     rep_rate_C = rep_rate*iter;
    
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
       

    %CREATION OF MATRIX FOR OUTPUT OF TABLE VALUES 
    fprintf('%3d...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f...%8.5f\n',...
        iter, theta1(1), theta1(2), theta1(3), theta1(4), theta1(5), theta1(6), theta1(7), theta1(8), theta1(9), theta1(10)); 

      
    repr = beta(tdata,theta1,m,a,b).*y(1,tdata,u,n,a,b)./((1-ExpCumDeathsN).*(gammasr + gammasd))...
          + (1 - alpha).*beta(tdata,theta1,m,a,b).*y(2,tdata,u,n,a,b)./((1-ExpCumDeathsN).*(gammavr + gammavd));     

    rp(:,iter) = repr;

end

tcalendar = datetime(2021,07,9) + caldays(1:length(tdata));

% legend([line10(1) line2a(1) line10x(1) line2x(1)],{'Reconstructed Transmission Rates','Mean of Reconstructed Transmission Rates',... 
%     'Effective Reproduction Number', 'Mean of Effective Reproduction Number'},'FontSize',18,'Location','best');
% 

figure
hold on
randColor = rand(NumCurves,3);
% cm = parula(NumCurves);
for it = 1:NumCurves
    txt=['Reporting Rate = ',num2str(rep_rate*it,'%.2f')];
    line10x = plot(tcalendar, rp(:,it),'Color',randColor(it,:),'DisplayName',txt);
% line10x = plot(tcalendar, rp(:,it),'Color',cm(it,:),'DisplayName',txt);
    set(line10x,'LineWidth',2)
    hold on
end
xlabel('\fontsize{18}Time (days)');
ylabel('\fontsize{18}R_e(t)')
ylim([0.0 2.1])
xlim([datetime(2021,07,9) datetime(2021,11,25)])

legend('FontSize',16,'Location','best') 
line5x = line([datetime(2021,07,9) datetime(2021,11,25)], [1 1],'Color',[0.6350, 0.0780, 0.1840]);
set(line5x,'LineWidth',2,'HandleVisibility', 'off')
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
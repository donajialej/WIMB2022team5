function curve_fitting_model1_paper

clear all
close all
clc

load AFluDat05-09.txt


tdata = AFluDat05_09(:,1);
qdata = AFluDat05_09(:,2);

tforward = (1:1:1825)';
format long

    function dy = model_1(t,y,k)


Lb = 1020/365;
L  = 1000/365;
mb = 1/(2*365);
nb = 0.1;
mu = 1/(65*365);
nu = 0.15; 
om = 127;

dy = zeros(5,1);   

dy(1) = Lb -(k(1)*sin(2*pi*(t+om)/365)+k(2))*y(1)*y(2)-mb*y(1);    
dy(2) = (k(1)*sin(2*pi*(t+om)/365)+k(2))*y(1)*y(2)-(nb+mb)*y(2);
dy(3) = L - k(3)*y(3)*y(2)-mu*y(3);
dy(4) = k(3)*y(3)*y(2)-(mu+nu)*y(4);
dy(5) = k(3)*y(3)*y(2);



    end

function q = model1(k,tdata)


[T,Y] = ode23s(@(t,y)(model_1(t,y,k)),tforward,[k(4)  k(5)  65000  k(6)  .00047]);

q = Y(tdata(:),5);

end

% 
% 
% %  k = [0.0000055  0.0000535  0.000000000023   1940.6   2.0698]; % Maia's new parameters
%  
% %  
% %  k= 1.0e+03 * [0.000000005104739   0.000000052574532   0.000000000000020...
% %     1.977392996478390   0.002084537225347];
% 
k =1.0e+03 *[0.000000005104570   0.000000052590476   0.000000000000019   1.976827190426774...
   0.002092395869505  0.000000351177667];
% k = 1.0e+03 *[0.000000005172004   0.000000052569586   0.000000000000020   1.977133326136062...
%               0.002086447001345   0.000000222735621]

%  k = 1.0e+03 *[0.000000003657819   0.000000050820140   0.000000000000025   2.040563930916215...
%      0.001294775681406   0.000000351177667] % best parameters from model 4
 
%  k = 1.0e+03 *[0.000000003630573   0.000000050787358   0.000000000000025   2.038636520883993...
%                0.001252000311376   0.000000017558883];

% k =1.0e+03 *[0.000000005050429   0.000000050798227   0.000000000000034   2.038636514548152...
%             0.001251999439234   0.000000017560865];
        
 lb = [0.0           0.0              0.0                   0.0      0.0       0.0];
%  
%  
% %  
  for i = 1:9
 
   [k,resnorm] = lsqcurvefit(@model1,k,tdata,qdata,lb,[],...
                             optimset('Disp','iter','TolX',10^(-15),'TolFun',10^(-15)))
  end

%   k =1.0e+03 *[0.000000005104570   0.000000052590509   0.000000000000019   1.976827190426774 0.002092395869505];

 figure1 = figure(2);
 
 axes1 = axes('Parent',figure1,...
    'YTickLabel',{'0.001','0.002','0.003','0.004','0.005'},...
    'YTick',[0.001 0.002 0.003 0.004 0.005],...
    'XTickLabel',{'2005','2006','2007','2008','2009','2010'},...
    'XTick',[0 365 730 1095 1460 1825],...
    'XGrid','on',...
    'AmbientLightColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
%     'GridLineStyle', '-',...
% set( hax , 'GridLineStyle' , '-' ) ;
% Uncomment the following line to preserve the X-limits of the axes
%  xlim(axes1,[0 1825]);
% Uncomment the following line to preserve the Y-limits of the axes
%  ylim(axes1,[0 0.005]);
box(axes1,'on');
hold(axes1,'all');

 [T,Y] = ode15s(@(t,y)(model_1(t,y,k)),tforward,[k(4)  k(5)  65000  .00047  .00047]);
 plot(tdata, qdata, 'Marker','.','Color',[1 0 0],...
                'MarkerSize',10,'LineStyle','none')
 plot(tforward, Y(:,5), '-b','LineWidth',1.2)
%  plot(tforward, Y(:,4), '-m','LineWidth',1.2)
 
 xlabel({'Time(days since 01/01/2005)'},'LineWidth',2,'FontSize',11,'FontName','Arial');
 ylabel({'Cumulative number of human cases (\times 10^5)'},'LineWidth',2,...
    'FontSize',11,...
    'FontName','Computer Modern');                      
                         

end


% Date 11/22

% k =
% 
%    1.0e+03 *
% 
%   Columns 1 through 4
% 
%    0.000000005172004   0.000000052569586   0.000000000000020   1.977133326136062
% 
%   Columns 5 through 6
% 
%    0.002086447001345   0.000000222735621
% 
% 
% resnorm =
% 
%      2.691310632784488e-07

% AIC  = 
% AICc = 

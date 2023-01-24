% Run after Param_Est_SVIRD_PCM_GA_normalized3.m
% The data can be loaded instead of re-runned

clear
close all
figures_for = 'DI1R'; % 'GA' 'CA' 'DI1R' 'DI2R' 
if strcmp(figures_for,'GA')
% Georgia:
savename_results = 'fits/phat_beta_data_estimates_GA.mat';
figurename_label = 'GA_';
elseif strcmp(figures_for,'CA')
% California:
savename_results = 'fits/phat_beta_data_estimates_CA.mat';
figurename_label = 'CA_';
elseif strcmp(figures_for,'DI1R')
% Synth DI1R
savename_results = 'fits/test/phat_beta_data_estimates_SynthDI1R.mat';
figurename_label = 'test/SynthDI1R_';
N = 10799566;
Y0N(1) = (N - 3942002 - 3580);
Y0N(2) = 3942002;
Y0N(3) = (3580 - 731);
Y0N(4) = 731;
Y0N(5) = 0;
Y0N(6) = 0;
figures_for = 'Synth';
elseif strcmp(figures_for,'DI2R')
% Synth DI1R
savename_results = 'fits/test/phat_beta_data_estimates_SynthDI2R.mat';
figurename_label = 'test/SynthDI2R_';
N = 39237836;
Y0N(1) = (N - 20086693 - 25039);
Y0N(2) = 20086693;
Y0N(3) = (25039 - 5116);
Y0N(4) = 5116;
Y0N(5) = 0;
Y0N(6) = 0;
figures_for = 'Synth';
end
load(savename_results)

if strcmp(figures_for,'Synth')
    [~,Y] = ode23s(@(t,y) svirdT(t,y,theta_model,N),tdata,Y0N);
end
%%    BETA RECONSTRUCTION FIGURE separate Reff graph

close all;
figure10 = figure;
colors_posible = [229,245,224; 161,217,155; 49,163,84; 166,206,227; 31,120,180]/255;
% betais, thetais, betamean, Reffis,Reffmean
t = tiledlayout(figure10,'flow');
for i = 1:size(Phatss,2)
ax1 = nexttile;
ax1.LineWidth = 0.8;
histogram(Phatss(:,i), 'FaceColor',	colors_posible(2,:),'EdgeAlpha',0.5)
ax1.FontSize = 6;
if i == 10
    xlabel('$\theta_{10}$','FontSize',8,'Interpreter','latex')
else
    xlabel('$\theta_'+string(i)+'$','FontSize',8,'Interpreter','latex')
end

if strcmp(figures_for,'Synth')
    linetheta = xline(theta_model(i),'--','LineWidth',1,'Color', [73,0,106]/255);
    xlim([(min(theta_model(i),min(Phatss(:,i)))-(max(Phatss(:,i))-min(Phatss(:,i)))/2) (max(theta_model(i),max(Phatss(:,i)))+(max(Phatss(:,i))-min(Phatss(:,i)))/2)])
end
   

end
%%% beta
ax1 = nexttile([1 5]);
ax1.LineWidth = 0.8;
hold on
line10 = plot(tcalendar, bt, 'Color',  colors_posible(1,:),'LineWidth',1,'HandleVisibility', 'off');
line2a = plot(tcalendar, mean(bt,2),'-', 'Color',  colors_posible(3,:),'LineWidth',1,'DisplayName','Mean $\beta(t)$');
if strcmp(figures_for,'Synth')
    plot(tcalendar, betamodel,'--','Color', [73,0,106]/255,'LineWidth',1,'DisplayName','True $\beta(t)$');
end
xlim([datetime(2021,07,9) datetime(2021,11,25)])
ax1.FontSize = 6;
leg = legend('FontSize',8,'Location','northeast', 'NumColumns',2,'Interpreter','latex');
leg.ItemTokenSize = [10,10];
%xlabel('Time (days)');
ylabel('$\beta(t)$','FontSize',8,'Interpreter','latex')
set(gca,'Xticklabel', []);


%%% Reff
ax1 = nexttile([1 5]);
ax1.LineWidth = 0.8;
hold on
line10x = plot(tcalendar, rp,'Color',  [colors_posible(4,:),0.5],'LineWidth',1,'HandleVisibility', 'off');
line5x = yline(1,':','LineWidth',1,'HandleVisibility', 'off');
line2x = plot(tcalendar,mean(rp,2),'-','Color',  colors_posible(5,:),'LineWidth',1,'DisplayName','Mean $$\mathcal{R}_e(t)$$');
if strcmp(figures_for,'Synth')
    plot(tcalendar, rpmodel,'--','Color',[73,0,106]/255,'LineWidth',1,'DisplayName','True $\mathcal{R}_e(t)$');
end
xlim([datetime(2021,07,9) datetime(2021,11,25)])
ax1.FontSize = 6;
leg = legend('FontSize',8,'Location','northeast', 'NumColumns',2,'Interpreter','latex');
leg.ItemTokenSize = [10,10];
%legend('boxoff')
xlabel('Time (days)','FontSize',8,'Interpreter','latex');
ylabel('$\mathcal{R}_e(t)$','FontSize',8,'Interpreter','latex');


t.Padding = 'none';
t.TileSpacing = 'none';
figure10.PaperUnits = 'centimeters';
figure10.PaperPosition = [.25 .25 11.7 10];
saveas(gcf,['figures/',figurename_label,'beta_reff'],'epsc')
savefig(['figures/',figurename_label,'beta_reff'])


%% Incidences and deaths

%%%%%% next figure
figure10 = figure;
t = tiledlayout(figure10,2,1);

ax1 = nexttile;
colors = [103,169,207; 1,108,89; 166,97,26]/255;
% curveis, mean curve, data 223,194,125 166,97,26
ax1.LineWidth = 0.8;
hold on
line1 = plot(tcalendar,curves_results,'Color',[colors(1,:),0.5],'LineWidth',0.8);
line4 = plot(tcalendar,mean(curves_results,2),'-','LineWidth',1,'Color',colors(2,:));
line3 = scatter(tcalendar,Idata, 50,'.','MarkerEdgeColor',colors(3,:));
ax1.FontSize = 6;
xlim([datetime(2021,07,9) datetime(2021,11,25)])
leg = legend([line4(1) line3(1)],{'Mean Incidence', ...
     [figures_for, ' Data']},'FontSize',8,'Location','northeast','Interpreter','latex');
leg.ItemTokenSize = [10,10];
ylabel('Case incidence','FontSize',8,'Interpreter','latex');
set(gca,'Xticklabel', []);

ax1 = nexttile;
colors = [117,107,177; 84,39,143; 230,97,1]/255;
% curveis, mean curve, data 
ax1.LineWidth = 0.8;
line1Z = plot(tcalendar,DC,'Color',[colors(1,:),0.3],'LineWidth',0.5);
hold on
line4Z = plot(tcalendar,mean(DC,2),'-','LineWidth',1,'Color',colors(2,:));
line3Z = scatter(tcalendar, Ddata, 50,'.','MarkerEdgeColor',colors(3,:));
ax1.FontSize = 6;
xlim([datetime(2021,07,9) datetime(2021,11,25)])
leg = legend([line4Z(1) line3Z(1) ],{'Mean $D(t)$',...
     [figures_for, ' Data']},'FontSize',8,'Location','northwest','Interpreter','latex');
leg.ItemTokenSize = [10,10];
xlabel('Time (days)','FontSize',8,'Interpreter','latex');
ylabel('Daily deaths','FontSize',8,'Interpreter','latex')

t.Padding = 'none';
t.TileSpacing = 'none';
figure10.PaperUnits = 'centimeters';
figure10.PaperPosition = [.25 .25 11.7 6];
saveas(gcf,['figures/',figurename_label,'data_fit'],'epsc')
savefig(['figures/',figurename_label,'data_fit'])


%% Recosntruction curves

figure10 = figure;
t = tiledlayout(figure10,2,1);

ax1 = nexttile;
colors = [29,145,192; 37,52,148; 1,108,89;1,70,54]/255;
% Sis, mean Sis, Vis, meanVis
ax1.LineWidth = 0.8;

line10A = plot(tcalendar, SC,'Color',[colors(1,:),0.5],'LineWidth',0.5,'HandleVisibility', 'off');
hold on
line2aA = plot(tcalendar,mean(SC,2),'-','LineWidth',1,'Color',colors(2,:),'DisplayName','Mean $S(t)$');
if strcmp(figures_for,'Synth')
    plot(tcalendar, Y(:,1),':','Color', [73,0,106]/255,'LineWidth',1,'DisplayName','True $S(t)$');
end
line10AB = plot(tcalendar, VC,'Color',[colors(3,:),0.5],'LineWidth',0.5,'HandleVisibility', 'off');
line2aAB = plot(tcalendar,mean(VC,2),'-','LineWidth',1,'Color',colors(4,:),'DisplayName','Mean $V(t)$');
if strcmp(figures_for,'Synth')
    plot(tcalendar, Y(:,2),'-.','Color', [73,0,106]/255,'LineWidth',1,'DisplayName','True $V(t)$');
end
ax1.FontSize = 6;

xlim([datetime(2021,07,9) datetime(2021,11,25)])
leg = legend('FontSize',8,'Location','northeast', 'NumColumns',2,'Interpreter','latex');
leg.ItemTokenSize = [10,10];
%xlabel('\fontsize{18}Time (days)');
ylabel('Num individuals','FontSize',8,'Interpreter','latex')
set(gca,'Xticklabel', []);

ax1 = nexttile;
colors = [203,24,29;165,15,21; 221,52,151; 174,1,126]/255;

% ISis, mean ISis, IVis, meanIVis
ax1.LineWidth = 0.8;
line10ABC = plot(tcalendar, ISC,'Color',[colors(1,:),0.5],'LineWidth',0.5,'HandleVisibility', 'off');
hold on
line2aABC = plot(tcalendar,mean(ISC,2),'-', 'LineWidth',1,'Color',colors(2,:),'DisplayName','Mean $I_s(t)$');
if strcmp(figures_for,'Synth')
    plot(tcalendar, Y(:,3),':','Color', [73,0,106]/255,'LineWidth',1,'DisplayName','True $I_s(t)$');
end
line10ABCD = plot(tcalendar, IVC,'Color',[colors(3,:),0.5],'LineWidth',0.5,'HandleVisibility', 'off');
line2aABCD = plot(tcalendar,mean(IVC,2),'-','LineWidth',1,'Color',colors(4,:),'DisplayName','Mean $I_v(t)$');
if strcmp(figures_for,'Synth')
    plot(tcalendar, Y(:,4),'-.','Color', [73,0,106]/255,'LineWidth',1,'DisplayName','True $I_v(t)$');
end
ax1.FontSize = 6;
%ylim([0 max(max(ISC))+10])
xlim([datetime(2021,07,9) datetime(2021,11,25)])


leg = legend('FontSize',8,'Location','northeast', 'NumColumns',2,'Interpreter','latex');
leg.ItemTokenSize = [10,10];
xlabel('Time (days)','FontSize',8,'Interpreter','latex');
ylabel('Num individuals','FontSize',8,'Interpreter','latex')


t.Padding = 'none';
t.TileSpacing = 'none';
figure10.PaperUnits = 'centimeters';
figure10.PaperPosition = [.25 .25 11.7 6];
saveas(gcf,['figures/',figurename_label,'reconstruction'],'epsc')
savefig(['figures/',figurename_label,'reconstruction'])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FUNCTION DEFINITIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = svirdT(t,y,theta,N)
m = 10;
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
dydt = zeros(6,1);
dydt(1) = -beta(t,theta,m,a,b).*y(1).*(y(3)+y(4))./(N - y(6)) - p.*y(1) + delta1*y(5) + delta2*y(2);
dydt(2) = p.*y(1) - (1-alpha)*beta(t,theta,m,a,b).*y(2).*(y(3)+y(4))./(N - y(6)) - delta2*y(2);
dydt(3) = beta(t,theta,m,a,b).*y(1).*(y(3)+y(4))./(N - y(6)) - (gammasr + gammasd)*y(3);
dydt(4) = (1-alpha).*beta(t,theta,m,a,b).*y(2).*(y(3)+y(4))./(N - y(6)) - (gammavr + gammavd)*y(4);
dydt(5) = gammasr.*y(3) + gammavr.*y(4) - delta1*y(5);
dydt(6) = gammasd.*y(3) + gammavd.*y(4);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yy = beta(t,theta,m,a,b) 
%global m a b 

       yy = 0;
   for j = 1:m
       yy = yy + theta(j).*legre(j-1,t,a,b);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = legre(j,t,A,B)

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

% Nice plots, run after modelversion1.m

close all;

colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]; ...
        [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]; [0.6350, 0.0780, 0.1840]; [0.8350, 0.6780, 0.1840]];
    
t = tiledlayout(2,1);
ax1 = nexttile;
for i =1:2
plot(T,y(:,i),'Color',colors(i,:),'LineWidth',2)
hold on
end
hold off
legend('S','V');
set(gca,'XTickLabel',[]);
ylabel('number')
ax2 = nexttile; 
for i =3:5
plot(T,y(:,i),'Color',colors(i,:),'LineWidth',2)
hold on
end
legend('I_s','I_v','D');
ylabel('number')
xlabel('days')


ax1.FontSize = 15;
ax2.FontSize = 15;
ax1.LineWidth = 1;
ax2.LineWidth = 1;
linkaxes([ax1,ax2],'x');
t.Padding = 'none';
t.TileSpacing = 'none';
saveas(gcf,'forward_model.png')

close all;
t = tiledlayout(2,1);
ax1 = nexttile;
plot(T,Inew,'Color',colors(2,:),'LineWidth',2)
hold on
plot(T,Inew_real,'*','Color',colors(1,:),'LineWidth',2)
hold off
legend('I_{new}','Noisy I_{new}');
set(gca,'XTickLabel',[]);
ylabel('number')
ax2 = nexttile; 
plot(T,y(:,5),'Color',colors(5,:),'LineWidth',2)
hold on
plot(T,D_real,'*','Color',colors(1,:),'LineWidth',2)
hold off
legend('D','Noisy D', 'Location','southeast');
ylabel('number')
xlabel('days')


ax1.FontSize = 15;
ax2.FontSize = 15;
ax1.LineWidth = 1;
ax2.LineWidth = 1;
linkaxes([ax1,ax2],'x');
t.Padding = 'none';
t.TileSpacing = 'none';
saveas(gcf,'synth_data.png')


close all;
ax1 = plot(T,beta,'Color',colors(1,:),'LineWidth',2);
legend('\beta(t)');
ylabel('\beta')
xlabel('days')
ylim([0,1])
set(gca,'Fontsize',20, 'LineWidth',1);
saveas(gcf,'beta_netherlands.png')


close all;
plot(T,Inew,'Color',colors(2,:),'LineWidth',2)
hold on
plot(T,Inew_real,'k*','LineWidth',2)
hold off
legend('I_{new}','Noisy I_{new}');
ylabel('number')
xlabel('days')
title('Netherlands data')
set(gca,'Fontsize',20, 'LineWidth',1);
saveas(gcf,'synth_data2.png')






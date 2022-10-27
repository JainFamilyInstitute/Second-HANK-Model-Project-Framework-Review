
mpar.maxlag=16;

signchange=1;


%% Dynare Plotting
folder='figure2/';

Name='baseline_3comp';
figurename=['IRF_Y_' com_fig_name '_' Name];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:1:16,IRF_Y(1:16),'b-','LineWidth',4.5)
hold on
plot(1:1:16,IRF_LIQUID_Y(1:16),'g--','LineWidth',4.5)
xlim([0 16])
ylim([-0.3 0.1])
set(gca,'XTick',[0 4 8 12 16])
set(gca,'YTick',[-0.3 -0.2 -0.1 0 0.1])
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot([0 16],[0 0],'--black') 
legend('HANK-2','HANK-1','Location','SouthEast')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.215 0.175 0.75 0.795])
printpdf(gcf, ['saves/' folder figurename])

figurename=['IRF_C_' com_fig_name '_' Name];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:1:16,IRF_C(1:16),'b-','LineWidth',4.5)
hold on
plot(1:1:16,IRF_LIQUID_C(1:16),'g--','LineWidth',4.5)
xlim([0 16])
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25) 
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
set(gca,'XTick',[0 4 8 12 16])
hold on
plot([0 16],[0 0],'--black') 
set(gca, 'FontName','arial','FontSize',25); 
set(gca,'Position',[0.25 0.175 0.7 0.795])
printpdf(gcf, ['saves/' folder figurename])

figurename=['IRF_I_' com_fig_name '_' Name];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:1:16,IRF_I(1:16),'b-','LineWidth',4.5)
hold on
plot(1:1:16,IRF_LIQUID_I(1:16),'g--','LineWidth',4.5)
xlim([0 16])
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25) 
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
set(gca,'XTick',[0 4 8 12 16])
hold on
plot([0 16],[0 0],'--black') 
set(gca, 'FontName','arial','FontSize',25); 
set(gca,'Position',[0.215 0.175 0.75 0.795])
printpdf(gcf, ['saves/' folder figurename])



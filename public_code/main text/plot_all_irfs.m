
mpar.maxlag=16;

signchange=1;


%% Dynare Plotting
folder='figure2/';

Name=['baseline_3comp'];
figurename=['IRF_Y_DynareNKIM_' Name];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot([1:1:16],IRF_Y(1:16),'b-','LineWidth',4.5)
hold on
plot([1:1:16],IRF_LIQUID_Y(1:16),'g--','LineWidth',4.5)
plot([1:1:16],signchange*100*y_eps_m(1:16),'r--o','LineWidth',4.5)
xlim([0 16])
ylim([-0.3 0.1])
set(gca,'XTick',[0 4 8 12 16])
set(gca,'YTick',[-0.3 -0.2 -0.1 0 0.1])
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot([0 16],[0 0],'--black') 
legend('HANK-2','HANK-1','RANK','Location','SouthEast')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.215 0.175 0.75 0.795])
printpdf(gcf, ['saves/' folder figurename])

figurename=['IRF_C_DynareNKIM_' Name];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot([1:1:16],IRF_C(1:16),'b-','LineWidth',4.5)
hold on
plot([1:1:16],IRF_LIQUID_C(1:16),'g--','LineWidth',4.5)
plot([1:1:16],signchange*100*c_eps_m(1:16),'r--o','LineWidth',4.5)
xlim([0 16])
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25) 
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
set(gca,'XTick',[0 4 8 12 16])
hold on
plot([0 16],[0 0],'--black') 
set(gca, 'FontName','arial','FontSize',25); 
set(gca,'Position',[0.25 0.175 0.7 0.795])
printpdf(gcf, ['saves/' folder figurename])

figurename=['IRF_I_DynareNKIM_' Name];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot([1:1:16],IRF_I(1:16),'b-','LineWidth',4.5)
hold on
plot([1:1:16],IRF_LIQUID_I(1:16),'g--','LineWidth',4.5)

plot([1:1:16],signchange*100*invest_eps_m(1:16),'r--o','LineWidth',4.5)
xlim([0 16])
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25) 
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
set(gca,'XTick',[0 4 8 12 16])
hold on
plot([0 16],[0 0],'--black') 
set(gca, 'FontName','arial','FontSize',25); 
set(gca,'Position',[0.215 0.175 0.75 0.795])
printpdf(gcf, ['saves/' folder figurename])


%%
folder='figure1/';


RR=(1/0.983-1)+1;
IRF_RANK_LP=100*100*((exp(q_eps_m(2:17))+((1+r_k_eps_m(2:17))*(RR-1)))./exp(q_eps_m(1:16))-RR)...
    -100*100*(R_eps_m(1:16)-pi_eps_m(2:17));

figurename=['IRF_RANK_ExpReturns_theta2_' num2str(100*par.theta_pi)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot([1:1:16],100*100*signchange*(R_eps_m(1:16)-pi_eps_m(2:17)),'r-*','LineWidth',4.5)
hold on
plot([1:1:16],100*100*((exp(q_eps_m(2:17))+((1+r_k_eps_m(2:17))*(RR-1)))./exp(q_eps_m(1:16))-RR),'g--','LineWidth',4.5)
plot(1:mpar.maxlag,zeros(1,16),'b-','LineWidth',4.5)

ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
set(gca,'XTick',[0 4 8 12 16])
ylim([-5 15])
xlim([0 16])

legend({'Liquid return','Illiquid return','Liquidity premium'},'Location','NorthEast','AutoUpdate','off')

plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.2125 0.175 0.75 0.795])
printpdf(gcf,['saves/' folder figurename])


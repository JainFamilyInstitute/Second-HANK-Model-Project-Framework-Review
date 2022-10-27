mpar.maxlag=16;
%% Plotting

figurename=['IRF_ExpReturns_theta2_' num2str(100*par.theta_pi)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,4*IRF_RBREAL(2:mpar.maxlag+1),'r-*','LineWidth',4.5)
hold on
plot(1:mpar.maxlag,4*IRF_RQ(1:mpar.maxlag),'g--','LineWidth',4.5)
plot(1:mpar.maxlag,4*IRF_LP(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
set(gca,'XTick',[0 4 8 12 16])
ylim([-20 60])
xlim([0 16])
legend({'Liquid return','Illiquid return','Liquidity premium'},'Location','NorthEast','AutoUpdate','off')
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.2125 0.175 0.75 0.795])
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_ReturnsPrices_theta2_' num2str(100*par.theta_pi)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,4*IRF_RBREAL(1:mpar.maxlag)./100,'b--','LineWidth',4.5)
hold on
plot(1:mpar.maxlag,4*IRF_RQ1(1:mpar.maxlag)./100,'r:','LineWidth',4.5)
plot(1:mpar.maxlag,IRF_W(1:mpar.maxlag),'g-.','LineWidth',4.5)
plot(1:mpar.maxlag,IRF_MARKUPS(1:mpar.maxlag),'k*','LineWidth',4.5)

ylabel('Percent or p.p.','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
set(gca,'XTick',[0 4 8 12 16])
% ylim([-5 15])
xlim([0 16])
legend({'Liquid return (p.p.)','Illiquid return (p.p.)','Wage (precent)','Profit (p.p.)'},'Location','SouthEast','AutoUpdate','off')
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.2125 0.175 0.75 0.795])
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_Y_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)   '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_Y(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
set(gca,'XTick',[0 4 8 12 16])
xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.2125 0.175 0.75 0.795])
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_C_theta2_' num2str(100*par.theta_pi)  '_rhoB_' num2str(100*par.rho_B)   '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_C(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_I_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_I(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.2 .2]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_G_theta2_' num2str(100*par.theta_pi)  '_rhoB_' num2str(100*par.rho_B)   '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_G(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.3 .1]);% set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_K_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_K(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_M_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_M(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-1 1]); 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_H_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)    '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_H(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_S_theta2_' num2str(100*par.theta_pi)  '_rhoB_' num2str(100*par.rho_B)   '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_S(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_RBPI_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_RB(1:mpar.maxlag),'r--','LineWidth',4.5)
hold on; 
plot(1:mpar.maxlag,IRF_RBREAL(1:mpar.maxlag),'b-','LineWidth',4.5)
hold on; 
legend({'nominal','real'},'Location','NorthEast','AutoUpdate','off')
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-1 1])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_RB_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_RB(2:mpar.maxlag+1),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-1 1]);
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_PI_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_PI(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
% %ylim([-1 1]);
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_Q_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_Q(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-1 1]);
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_D_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_D(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-1 1]);
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_W_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_W(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_N_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_N(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_PF_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_PROFITS(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
%ylim([-.1 .1]); set(gca,'YTick',[-.1 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_LP_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag,IRF_LP(1:mpar.maxlag),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
% set(gca,'YTick',[-.2 -.1 0]);
% set(gca,'YTick',[-.1 0 .1]) 
%ylim([-.1 0.1]);
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])

figurename=['IRF_MK_theta2_' num2str(100*par.theta_pi) '_rhoB_' num2str(100*par.rho_B)     '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
IRF_MK=100*log((M./K)/(grid.B/grid.K));%(100*(M./K-grid.B/grid.K));
plot(1:mpar.maxlag,IRF_MK(1:mpar.maxlag),'r-o','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
% set(gca,'YTick',[-.2 -.1 0]);
% set(gca,'YTick',[-.1 0 .1]) 
%ylim([-.1 0.1]);
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')

set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['../saves/' figurename])



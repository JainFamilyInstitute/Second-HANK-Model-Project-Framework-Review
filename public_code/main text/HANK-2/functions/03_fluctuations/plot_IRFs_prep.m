
IRF_distr=Gamma_state*IRF_state_sparse(1:numstates-2,1:mpar.maxlag);
% preparation
IRF_H=100*grid.h(1:end-1)*IRF_distr(mpar.nm+mpar.nk+(1:mpar.nh-1),2:end)/par.H;
K=grid.k*IRF_distr(mpar.nm+(1:mpar.nk),:)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*grid.k*IRF_distr(mpar.nm+(1:mpar.nk),2:end)./grid.K;
IRF_M=100*grid.m*IRF_distr((1:mpar.nm),2:end)./(targets.B+par.ABS*grid.K);
M=grid.m*IRF_distr((1:mpar.nm),1:end)+targets.B-par.ABS*(K(1:end)-grid.K);
IRF_S=100*IRF_state_sparse(n2(2)+2,1:end-1);
IRF_Y=100*IRF_state_sparse(end-oc+3,1:end-1);
IRF_G=100*IRF_state_sparse(end-oc+4,1:end-1);
IRF_T=100*IRF_state_sparse(end-oc+9,1:end-1);

IRF_W=100*IRF_state_sparse(end-oc+5,1:end-1);
W=exp(log(par.W)+IRF_state_sparse(end-oc+5,1:end-1));
IRF_N=100*IRF_state_sparse(end-oc+8,1:end-1);
N=par.N*(1+IRF_state_sparse(end-oc+8,1:end-1));
IRF_PROFITS=100*IRF_state_sparse(end-oc+7,1:end-1);
PROFITS=exp(log(par.PROFITS)+IRF_state_sparse(end-oc+7,1:end-1));

MARKUPS=exp(log(1-par.mu)+IRF_state_sparse(end-oc+13,1:end-1));
IRF_MARKUPS=100*(MARKUPS-(1-par.mu));
IRF_C=100*IRF_state_sparse(end-oc+12,1:end-1);

IRF_R=100*IRF_state_sparse(end-oc+6,1:end-1);
IRF_PI=100*100*IRF_state_sparse(end-oc+2,1:end-1);
PI=1+IRF_state_sparse(end-oc+2,1:end-1);
Q=exp(log(par.Q)+IRF_state_sparse(end-oc+1,1:end-1));
R=exp(log(par.R)+IRF_state_sparse(end-oc+6,1:end-1));
RB=exp(log(par.RB)+(IRF_state_sparse(n2(2)+1,1:end-1)));
IRF_RB=100*100*(RB-par.RB);

IRF_RBREAL=100*100*(RB(1:end)./PI(1:end)-par.RB);

IRF_Q=100*100*(Q-par.Q);
IRF_D=100*100*((1+IRF_R/100)*par.R-par.R);

IRF_RQ=100*100*((Q(2:end)+R(2:end))./Q(1:end-1)-(1+par.R));
IRF_RQ1=100*100*((Q(1:end)+R(1:end))./[1 Q(1:end-1)]-(1+par.R));

IRF_LP=100*100*(((Q(2:end)+R(2:end))./[Q(1:end-1)]-RB(2:end)./PI(2:end))-((1+par.R/par.Q)-par.RB));
MK=grid.B./grid.K;

IRFname=['../saves/IRF_' casename '_theta2_' num2str(100*par.theta_pi) '_theta3_' num2str(100*par.rho_R)  '_gamma2_' num2str(100*par.gamma_T) '_ABS_' num2str(par.ABS*100) '_' aggrshock];
save(IRFname,'IRF_Y','IRF_C','IRF_I','IRF_G','IRF_M','IRF_K','IRF_Q','IRF_N','IRF_D','IRF_PI','IRF_W','IRF_PROFITS','IRF_LP','IRF_RB','IRF_S','K','M','MK');

%% Plotting
mpar.maxlag=16;

folder='figure1/';

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
printpdf(gcf,['../saves/' folder figurename])





folder='figure3/';


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
printpdf(gcf,['../saves/' folder figurename])


close all

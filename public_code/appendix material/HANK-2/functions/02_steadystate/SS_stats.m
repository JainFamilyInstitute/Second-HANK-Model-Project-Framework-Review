% Steady state statistics
%%

[mesh.k,mesh.m]=meshgrid(grid.k,grid.m);
clear targets
targets.ShareBorrower=sum((grid.m<0).*sum(sum(joint_distr,2),3)');
targets.K=sum(grid.k.*sum(sum(joint_distr,1),3));
targets.B=grid.m*sum(sum(joint_distr,2),3);
grid.K=targets.K;
grid.B=targets.B;

JDredux = sum(joint_distr,3);
targets.BoverK     = targets.B/targets.K;

targets.L=grid.N*sum(grid.h'.*squeeze(sum(sum(joint_distr,1),2)));
targets.KY=targets.K/Output;
targets.BY=targets.B/Output;
targets.Y=Output;
BCaux_M=sum(sum(joint_distr,2),3);
targets.m_bc=BCaux_M(1,:);
targets.m_0=BCaux_M(grid.m==0);
BCaux_K=sum(sum(joint_distr,1),3);
targets.k_bc=BCaux_K(:,1);
aux_MK=sum(joint_distr,3);
targets.WtH_b0=sum(aux_MK(mesh.m==0 & mesh.k>0));
targets.WtH_bnonpos=sum(aux_MK(mesh.m<=0 & mesh.k>0));

TN=(1-par.tau).*(n_n_guess).*W_fc.*meshes.h./par.H;
TA=(1-par.tau).*(n_a_guess).*W_fc.*meshes.h./par.H;
TT=[TA(:); TN(:)];
TT=[par.nu.*joint_distr(:);(1-par.nu).*joint_distr(:)]'*TT(:);
targets.T =TT +(1-par.tau).*Profits_fc;
par.G=targets.B*(1-par.RB/par.PI)+targets.T;
par.R=R_fc;
par.W=W_fc(1);
par.PROFITS=Profits_fc(1);
par.N=grid.N;
targets.GtoY=par.G/Output;

c=[par.nu.*joint_distr(:);(1-par.nu).*joint_distr(:)]'*[c_a_guess(:); c_n_guess(:)];
targets.C = c;
par.spread=par.R-par.RB+1;

%% Ginis
% Net worth Gini
mplusk=mesh.k(:)*par.Q+mesh.m(:);
[mplusk, IX]       = sort(mplusk);
moneycapital_pdf   = JDredux(IX);
moneycapital_cdf   = cumsum(moneycapital_pdf);
targets.NegNetWorth= sum((mplusk<0).*moneycapital_pdf);

S                  = cumsum(moneycapital_pdf.*mplusk)';
S                  = [0 S];
targets.GiniW      = 1-(sum(moneycapital_pdf.*(S(1:end-1)+S(2:end))')/S(end));

% Liquid Gini
[liquid_sort, IX]  = sort(mesh.m(:));
liquid_pdf         = JDredux(IX);
liquid_cdf         = cumsum(liquid_pdf);
targets.Negliquid  = sum((liquid_sort<0).*liquid_pdf);

S                  = cumsum(liquid_pdf.*liquid_sort)';
S                  = [0 S];
targets.GiniLI      = 1-(sum(liquid_pdf.*(S(1:end-1)+S(2:end))')/S(end));

% Illiquid Gini
[illiquid_sort, IX] = sort(mesh.k(:));
illiquid_pdf        = JDredux(IX);
illiquid_cdf        = cumsum(illiquid_pdf);
targets.Negliquid   = sum((illiquid_sort<0).*illiquid_pdf);

S                   = cumsum(illiquid_pdf.*illiquid_sort)';
S                   = [0 S];
targets.GiniIL      = 1-(sum(illiquid_pdf.*(S(1:end-1)+S(2:end))')/S(end));



%% Liquidity
folder='figure2\';

c_0=par.nu.*c_a_guess+(1-par.nu).*c_n_guess;
m_0=par.nu.*m_a_star+(1-par.nu).*m_n_star;
k_0=par.nu.*cap_a_star+(1-par.nu).*meshes.k;

wealthsort='LiquidWealth';

[lambdaLIP,lambdaILP,lambdaC,lambdaLI,lambdaIL,locstdIL,locstdLI,lambdaLIINC] = diststats_liquidity(joint_distr,c_0,m_0,k_0,grid,mpar,par,par.Q,wealthsort);

figurename=['Model_Average_LiquidityBAR_' casename '_by_' wealthsort];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
aux=(lambdaLI(1:end,:)./lambdaIL(1:end,:))*100;
LIRATIOSQUINTILES=mean(reshape(aux,[20 5]),1);
bar(1:1:5,LIRATIOSQUINTILES);
ylabel('Liquid to illiquid assets (%)','Interpreter','none','FontName','arial','FontSize',20)
xlabel('Quintile of liquid wealth','Interpreter','none','FontName','arial','FontSize',20)
hold on
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.125 0.165 0.78 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

wealthsort='NetWealth';
[lambdaLIP,lambdaILP,lambdaC,lambdaLI,lambdaIL,locstdIL,locstdLI,lambdaLIINC] = diststats_liquidity(joint_distr,c_0,m_0,k_0,grid,mpar,par,par.Q,wealthsort);

figurename=['Model_Average_LiquidityBAR_' casename '_by_' wealthsort];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
aux=(lambdaLI(1:end,:)./lambdaIL(1:end,:))*100;
LIRATIOSQUINTILES=mean(reshape(aux,[20 5]),1);
bar(1:1:5,LIRATIOSQUINTILES);
ylabel('Liquid to illiquid assets (%)','Interpreter','none','FontName','arial','FontSize',20)
xlabel('Quintile of liquid wealth','Interpreter','none','FontName','arial','FontSize',20)
hold on
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.125 0.165 0.78 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

targets.LIratiomedian=mean(lambdaLI([51],:)./lambdaIL([51],:));

targets

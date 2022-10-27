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

%% Produce Figure 4
wealthsort='LiquidWealth';
folder='figure4\';
MPS_a_m = zeros(mpar.nm,mpar.nk,mpar.nh);
MPS_n_m = zeros(mpar.nm,mpar.nk,mpar.nh);
for kk=1:mpar.nk
    for hh=1:mpar.nh
        MPS_a_m(:,kk,hh)=gradient(squeeze(m_a_star(:,kk,hh)))./gradient(grid.m)';
        MPS_n_m(:,kk,hh)=gradient(squeeze(m_n_star(:,kk,hh)))./gradient(grid.m)';
    end
end

[lambdaMPS_a_m,lambdaMPS_n_m] = diststats_liquidity(joint_distr,MPS_a_m,MPS_a_m,MPS_n_m,grid,mpar,par,1,wealthsort);

% figurename=['Model_MPS_' casename];
% figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
% plot(1:100,100*(1-lambdaMPS_a_m),'b-','LineWidth',4.5)
% hold on
% plot(1:100,100*(1-lambdaMPS_n_m),'r--','LineWidth',4.5)
% legend({'adjusters','non-adjusters'},'Location','NorthWest')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
% xlabel('Percentile of liquid wealth','Interpreter','none','FontName','arial','FontSize',25)
% set(gca, 'FontName','arial','FontSize',20);
% set(gca,'Position',[0.15 0.175 0.8 0.8])
% printpdf(gcf, ['..\saves\' folder figurename])


%% Produce Figure 4
MPC_a_m = zeros(mpar.nm,mpar.nk,mpar.nh);
MPC_n_m = zeros(mpar.nm,mpar.nk,mpar.nh);
for kk=1:mpar.nk
    for hh=1:mpar.nh
        MPC_a_m(:,kk,hh)=gradient(squeeze(c_a_guess(:,kk,hh)))./gradient(grid.m)';
        MPC_n_m(:,kk,hh)=gradient(squeeze(c_n_guess(:,kk,hh)))./gradient(grid.m)';
    end
end


[lambdaMPC_a_m,lambdaMPC_n_m] = diststats_liquidity(joint_distr,MPC_a_m,MPC_a_m,MPC_n_m,grid,mpar,par,1,wealthsort);

figurename=['Model_MPC_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(1:100,100*lambdaMPC_a_m,'b-','LineWidth',4.5)
hold on
plot(1:100,100*lambdaMPC_n_m,'r--','LineWidth',4.5)
legend({'adjusters','non-adjusters'},'Location','NorthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of liquid wealth','Interpreter','none','FontName','arial','FontSize',25)
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.15 0.175 0.8 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

%% Produce Figure 4
MPI_a_m = zeros(mpar.nm,mpar.nk,mpar.nh);
for kk=1:mpar.nk
    for hh=1:mpar.nh
        MPI_a_m(:,kk,hh)=gradient(squeeze(cap_a_star(:,kk,hh)))./gradient(grid.m)';
    end
end

[lambdaMPI_a_m] = diststats_liquidity(joint_distr,MPI_a_m,MPI_a_m,MPI_a_m,grid,mpar,par,1,wealthsort);

figurename=['Model_MPI_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(1:100,100*lambdaMPI_a_m,'b-','LineWidth',4.5)
hold on
plot(1:100,zeros(1,100),'r--','LineWidth',4.5)
legend({'adjusters','non-adjusters'},'Location','SouthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of liquid wealth','Interpreter','none','FontName','arial','FontSize',25)
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.15 0.175 0.8 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

%% Produce Figure 4
figurename=['Model_MPCplusMPI_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(1:100,100*(lambdaMPI_a_m+lambdaMPC_a_m),'b-','LineWidth',4.5)
hold on
plot(1:100,100*(lambdaMPC_n_m),'r--','LineWidth',4.5)
legend({'adjusters','non-adjusters'},'Location','NorthWest')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of liquid wealth','Interpreter','none','FontName','arial','FontSize',25)
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.15 0.175 0.8 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

figurename=['Model_MPCMPIMPSall_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(1:100,100*(par.nu*lambdaMPC_a_m+(1-par.nu)*lambdaMPC_n_m),'b-','LineWidth',4.5)
hold on
plot(1:100,100*(par.nu*lambdaMPI_a_m+(1-par.nu)*0),'r--','LineWidth',4.5)
plot(1:100,100*(par.nu*((lambdaMPI_a_m+lambdaMPC_a_m))+(1-par.nu)*(lambdaMPC_n_m)),'g-.','LineWidth',4.5)
legend({'MPC','MPI','MPS'},'Location','NorthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of liquid wealth','Interpreter','none','FontName','arial','FontSize',25)
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.15 0.175 0.8 0.8])
printpdf(gcf, ['..\saves\' folder figurename])


targets

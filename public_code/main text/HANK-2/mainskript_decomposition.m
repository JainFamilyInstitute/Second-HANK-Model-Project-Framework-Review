
%% Produce IRFs
mpar.maxlag=1000;

x0=zeros(numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

PI=1+IRF_state_sparse(end-oc+2,1:end-1);
Q=exp(log(par.Q)+IRF_state_sparse(end-oc+1,1:end-1));
R=exp(log(par.R)+IRF_state_sparse(end-oc+6,1:end-1));
RB=exp(log(par.RB)+(IRF_state_sparse(n2(2)+1,1:end-1)));
W=exp(log(par.W)+IRF_state_sparse(end-oc+5,1:end-1));
PROFITS=exp(log(par.PROFITS)+IRF_state_sparse(end-oc+7,1:end-1));

%% MIT Shock
mutil = @(c)(1./(c.^par.xi));
[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

meshaux=meshes;
meshaux.h(:,:,end)=1000;
TT=mpar.maxlag;
casename            = 'SS_BASELINE_ALL'

IRF_MIT=zeros(TT-1,4,2);

for j=1:4
    
rT=ones(TT)*par.R;
QT=ones(TT);
RBT=ones(TT)*par.RB;
PIT=ones(TT);
WT=ones(TT)*par.W;
PROT=ones(TT)*par.PROFITS;

if j==1
    RBT=[par.RB RB];
    PIT=[1 PI];
elseif j==2
    WT=[par.W W];
elseif j==3
    PROT=[par.PROFITS PROFITS];
elseif j==4
    rT=[par.R R];
    QT=[1 Q];
end

ZT=ones(size(rT));
KT = ones(size(rT))*targets.K;
NT = ones(size(rT))*par.N;

MPAT(:,:,:,1) = m_a_star;
MPNT(:,:,:,1) = m_n_star;
CPAT(:,:,:,1) = c_a_guess;
CPNT(:,:,:,1) = c_n_guess;
VkT(:,:,:,1) = Vk;
KPT(:,:,:,1) = cap_a_star;

MPAT(:,:,:,TT) = m_a_star;
MPNT(:,:,:,TT) = m_n_star;
CPAT(:,:,:,TT) = c_a_guess;
CPNT(:,:,:,TT) = c_n_guess;
VkT(:,:,:,TT) = Vk;
KPT(:,:,:,TT) = cap_a_star;

% Backwards iteration of policy
for h = 1:(TT-2)
    t = TT-h;
 
    % Incomes (grids)
    inc.rent    = meshes.k*rT(t);
    inc.capital = meshes.k*QT(t);
    inc.money   = (RBT(t)/PIT(t)).*meshes.m...
        + (meshes.m<0).*(par.borrwedge/PIT(t)).*meshes.m;

    mutil_c =  par.nu*mutil(CPAT(:,:,:,t+1))+(1-par.nu)*mutil(CPNT(:,:,:,t+1));
    Vk = VkT(:,:,:,t+1);

  
    % First Set: Value Functions, Marginal Values
    % Update policies    
    EVk = reshape(reshape(Vk,[mpar.nm*mpar.nk mpar.nh])*P_H',[mpar.nm,mpar.nk mpar.nh]);
    RBaux = RBT(t+1)/PIT(t+1) + (meshes.m<0).*(par.borrwedge/PIT(t+1));   
    EVm = reshape(reshape(RBaux(:).*mutil_c(:),[mpar.nm*mpar.nk mpar.nh])*P_H',[mpar.nm,mpar.nk mpar.nh]);
    
    [c_a_t,m_a_t,k_t,c_n_t,m_n_t,nnt,nat] = EGM_policyupdate_MIT(EVm,EVk,QT(t),PIT(t),RBT(t),WT(t),PROT(t),par.tau,inc,meshes,grid,par,mpar);

    naminus= ((1/par.gamma).* mutil(c_a_t(:)).*meshes.h(:).*par.tau.*WT(t)).^(par.sigma_n);
    naminus          = min(naminus,1);
    naminus          = reshape(naminus,[mpar.nm mpar.nk mpar.nh]);
    naminus(:,:,end) = 0;
    nnminus= ((1/par.gamma).* mutil(c_n_t(:)).*meshes.h(:).*par.tau.*WT(t)).^(par.sigma_n);
    nnminus          = min(nnminus,1);
    nnminus          = reshape(nnminus,[mpar.nm mpar.nk mpar.nh]);
    nnminus(:,:,end) = 0;
  
    Vk_next = griddedInterpolant(meshaux.m,meshaux.k,meshaux.h,reshape(EVk,[mpar.nm mpar.nk mpar.nh]));
    Vk_aux  = par.nu.*(rT(t)+QT(t)).*mutil(c_a_t) + (1-par.nu).*rT(t).*mutil(c_n_t)+ par.beta.*(1-par.nu).*Vk_next(m_n_t,meshaux.k,meshaux.h); % Expected marginal utility at consumption policy (w &w/o adjustment)
    VkT(:,:,:,t)=Vk_aux;
    
    CPAT(:,:,:,t)=c_a_t;    
    CPNT(:,:,:,t)=c_n_t;
    KPT(:,:,:,t) = k_t;
    MPAT(:,:,:,t) = m_a_t;
    MPNT(:,:,:,t) = m_n_t;
   
   
end

POLT(:,:,:,:,1,j)=CPAT(:,:,:,1:2);
POLT(:,:,:,:,2,j)=CPNT(:,:,:,1:2);
POLT(:,:,:,:,3,j)=MPAT(:,:,:,1:2);
POLT(:,:,:,:,4,j)=MPNT(:,:,:,1:2);
POLT(:,:,:,:,5,j)=KPT(:,:,:,1:2);

% Forward Iteration of Distribution
DIST(:,1)  = joint_distr(:);
DIST(:,TT) = joint_distr(:);

for t=2:TT

    [H]=Gen_BigTransH(MPAT(:,:,:,t-1), MPNT(:,:,:,t-1),KPT(:,:,:,t-1), P_H,par, mpar, grid);
    
    DIST(:,t)=DIST(:,t-1)'*H;

end

% PARTAL IRFS
K_MIT=zeros(TT,1);
C_MIT=zeros(TT,1);

for t=1:TT

    
    KAaux=KPT(:,:,:,t);
    KNaux=meshes.k;    
    
    CAaux=CPAT(:,:,:,t);
    CNaux=CPNT(:,:,:,t);
    Daux=DIST(:,t);
    
    K_MIT(t)=par.nu*Daux(:)'*KAaux(:)+(1-par.nu).*Daux(:)'*KNaux(:);

    C_MIT(t)=par.nu*Daux(:)'*CAaux(:)+(1-par.nu).*Daux(:)'*CNaux(:);

end

Iaux=(K_MIT(2:end)-(1-par.delta)*K_MIT(1:end-1));
I_MIT=100*(Iaux/(par.delta*grid.K)-1);
C_MIT_IRF=100*(C_MIT/C_MIT(1)-1);

IRF_MIT(:,j,1)=I_MIT;
IRF_MIT(:,j,2)=C_MIT_IRF(2:end);

end


%% Produce Figure 3
mpar.maxlag=18;

IRF_I_MIT=IRF_MIT(1:mpar.maxlag-2,:,1);
IRF_C_MIT=IRF_MIT(1:mpar.maxlag-2,:,2);

folder='figure3\';
figurename=['IRF_I_DECOMP_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag-2,IRF_I(1:mpar.maxlag-2),'k-','LineWidth',4.5)
hold on
plot(1:mpar.maxlag-2,IRF_I_MIT(:,1),'b--','LineWidth',4.5)
plot(1:mpar.maxlag-2,IRF_I_MIT(:,4),'r:','LineWidth',4.5)
plot(1:mpar.maxlag-2,IRF_I_MIT(:,2),'g-.','LineWidth',4.5)
plot(1:mpar.maxlag-2,IRF_I_MIT(:,3),'k*','LineWidth',4.5)

legend({'Total','Liquid return','Illiquid return','Wage','Profit'},'Location','SouthEast','AutoUpdate','off')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['..\saves\' folder figurename])


figurename=['IRF_C_DECOMP_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
plot(1:mpar.maxlag-2,IRF_C(1:mpar.maxlag-2),'k-','LineWidth',4.5)
hold on
plot(1:mpar.maxlag-2,IRF_C_MIT(:,1),'b--','LineWidth',4.5)
plot(1:mpar.maxlag-2,IRF_C_MIT(:,4),'r:','LineWidth',4.5)
plot(1:mpar.maxlag-2,IRF_C_MIT(:,2),'g-.','LineWidth',4.5)
plot(1:mpar.maxlag-2,IRF_C_MIT(:,3),'k*','LineWidth',4.5)

legend({'Total','Liquid return','Illiquid return','Wage','Profit'},'Location','SouthEast','AutoUpdate','off')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on; set(gca,'XTick',[0 4 8 12 16]); xlim([0 16])
ylim([-.2 .1]); set(gca,'YTick',[-.2 -.1 0 .1]) 
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2125 0.175 0.75 0.795]) 
printpdf(gcf,['..\saves\' folder figurename])



%% Decomposition across wealth distribution
for j=1:4
    
x_A0=POLT(:,:,:,1,1,j);%CPAT(:,:,:,1);
m_A0=POLT(:,:,:,1,3,j);%MPAT(:,:,:,1);
k_A0=POLT(:,:,:,1,5,j);%KPT(:,:,:,1);

x_N0= POLT(:,:,:,1,2,j);% CPNT(:,:,:,1);
m_N0=POLT(:,:,:,1,4,j);%MPNT(:,:,:,1);
k_N0=meshes.k;

[lambdaLI_A0(:,j),lambdaIL_A0(:,j),lambdaC_A0(:,j)] = diststats_liquidity_weights(joint_distr(:),m_A0,k_A0,x_A0,grid,1,1,mpar);
[lambdaLI_N0(:,j),lambdaIL_N0(:,j),lambdaC_N0(:,j)] = diststats_liquidity_weights(joint_distr(:),m_N0,k_N0,x_N0,grid,1,1,mpar);

lambdaLI_0(:,j)=par.nu * lambdaLI_A0(:,j) + (1-par.nu) * lambdaLI_N0(:,j);
lambdaIL_0(:,j)=par.nu * lambdaIL_A0(:,j) + (1-par.nu) * lambdaIL_N0(:,j);
lambdaC_0(:,j)=par.nu * lambdaC_A0(:,j) + (1-par.nu) * lambdaC_N0(:,j);

x_A1=POLT(:,:,:,2,1,j);%CPAT(:,:,:,1);
m_A1=POLT(:,:,:,2,3,j);%MPAT(:,:,:,1);
k_A1=POLT(:,:,:,2,5,j);%KPT(:,:,:,1);

x_N1=POLT(:,:,:,2,2,j);% CPNT(:,:,:,1);
m_N1=POLT(:,:,:,2,4,j);%MPNT(:,:,:,1);
k_N1=meshes.k;

[lambdaLI_A1(:,j),lambdaIL_A1(:,j),lambdaC_A1(:,j)] = diststats_liquidity_weights(joint_distr(:),m_A1,k_A1,x_A1,grid,1,1,mpar);
[lambdaLI_N1(:,j),lambdaIL_N1(:,j),lambdaC_N1(:,j)] = diststats_liquidity_weights(joint_distr(:),m_N1,k_N1,x_N1,grid,1,1,mpar);

lambdaLI_1(:,j)=par.nu * lambdaLI_A1(:,j) + (1-par.nu) * lambdaLI_N1(:,j);
lambdaIL_1(:,j)=par.nu * lambdaIL_A1(:,j) + (1-par.nu) * lambdaIL_N1(:,j);
lambdaC_1(:,j)=par.nu * lambdaC_A1(:,j) + (1-par.nu) * lambdaC_N1(:,j);

end
%% Produce Figure 5
cutp=11;
folder='figure5\';


auxirf=100*(lambdaC_1(1:end,:)./lambdaC_0(1:end,:)-1);
figurename=['IRF_ConsumptionMIT_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(1:100,sum(auxirf,2),'k-','LineWidth',4)
hold on
plot(1:100,auxirf(1:end,1),'b--','LineWidth',4)
plot(1:100,auxirf(1:end,4),'r:','LineWidth',4)
plot(1:100,auxirf(1:end,2),'g-.','LineWidth',4)
plot(1:100,auxirf(1:end,3),'k*','LineWidth',0.25)
legend({'Total','Liquid return','Illiquid return','Wage','Profit'},'Location','SouthEast','AutoUpdate','off')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',20)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',20)
hold on
xlim([1 100])

plot(1:100,zeros(1,100),'--black')
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.175 0.175 0.78 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

auxirf=100*(lambdaLI_1(cutp:end,:)./lambdaLI_0(cutp:end,:)-1);
figurename=['IRF_MoneyMIT_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(cutp:100,sum(auxirf,2),'k-','LineWidth',4)
hold on
plot(cutp:100,auxirf(1:end,1),'b--','LineWidth',4)
plot(cutp:100,auxirf(1:end,4),'r:','LineWidth',4)
plot(cutp:100,auxirf(1:end,2),'g-.','LineWidth',4)
plot(cutp:100,auxirf(1:end,3),'k*','LineWidth',0.25)
legend({'Total','Liquid return','Illiquid return','Wage','Profit'},'Location','SouthEast','AutoUpdate','off')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',20)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',20)
hold on
xlim([cutp 100])

plot(cutp:100,zeros(1,100-cutp+1),'--black')
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.175 0.175 0.78 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

auxirf=100*(lambdaIL_1(cutp:end,:)./lambdaIL_0(cutp:end,:)-1);
figurename=['IRF_CapitalMIT_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(cutp:100,sum(auxirf,2),'k-','LineWidth',4)
hold on
plot(cutp:100,auxirf(1:end,1),'b--','LineWidth',4)
plot(cutp:100,auxirf(1:end,4),'r:','LineWidth',4)
plot(cutp:100,auxirf(1:end,2),'g-.','LineWidth',4)
plot(cutp:100,auxirf(1:end,3),'k*','LineWidth',0.25)
legend({'Total','Liquid return','Illiquid return','Wage','Profit'},'Location','SouthEast','AutoUpdate','off')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',20)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',20)
hold on
xlim([cutp 100])

plot(cutp:100,zeros(1,100-cutp+1),'--black')
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.175 0.175 0.78 0.8])
printpdf(gcf, ['..\saves\' folder figurename])

auxirf=100*((lambdaLI_1(cutp:end,:)./lambdaIL_1(cutp:end,:))-(lambdaLI_0(cutp:end,:)./lambdaIL_0(cutp:end,:)));
figurename=['IRF_PortLiqPPTMIT_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 500])
plot(cutp:100,sum(auxirf,2),'k-','LineWidth',4)
hold on
plot(cutp:100,auxirf(1:end,1),'b--','LineWidth',4)
plot(cutp:100,auxirf(1:end,4),'r:','LineWidth',4)
plot(cutp:100,auxirf(1:end,2),'g-.','LineWidth',4)
plot(cutp:100,auxirf(1:end,3),'k*','LineWidth',0.25)
legend({'Total','Liquid return','Illiquid return','Wage','Profit'},'Location','SouthEast','AutoUpdate','off')
ylabel('Percentage points','Interpreter','none','FontName','arial','FontSize',20)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',20)
hold on
xlim([cutp 100])
plot(cutp:100,zeros(1,100-cutp+1),'--black')
set(gca, 'FontName','arial','FontSize',20);
set(gca,'Position',[0.175 0.175 0.78 0.8])
printpdf(gcf, ['..\saves\' folder figurename])



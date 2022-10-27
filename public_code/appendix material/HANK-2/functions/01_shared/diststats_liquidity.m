function [lambdaLIP,lambdaILP,lambdaC,lambdaLI,lambdaIL,locstdIL,locstdLI,lambdaLIINC] = diststats_liquidity(joint_distr,c_prime,m_prime,k_prime,grid,mpar,par,Q,wealthsort)
% function [quints,gini,constrained] = diststats(joint_distr,kgrid)
% distats computes the following summary statistics given 
% distribution of agents over productivity and assets:
% - asset quintiles
% - gini coefficient
% - share of borrowing constrained agents
%
% author: Ralph Lutticke
% Date of last change: May 2018

joint_distr=reshape(joint_distr, [mpar.nm mpar.nk mpar.nh]);

%% Networth
[mesh3.k,mesh3.m,mesh3.h]=meshgrid(grid.k,grid.m,grid.h);
for j=1:mpar.nm
    for i=1:mpar.nk
        Inc_condh(j,i)=squeeze(joint_distr(j,i,:)./sum(joint_distr(j,i,:)))' * squeeze(mesh3.h(j,i,:));
        m_policy(j,i)=squeeze(joint_distr(j,i,:)./sum(joint_distr(j,i,:)))' * squeeze(m_prime(j,i,:));
        k_policy(j,i)=squeeze(joint_distr(j,i,:)./sum(joint_distr(j,i,:)))' * squeeze(k_prime(j,i,:));
        c_policy(j,i)=squeeze(joint_distr(j,i,:)./sum(joint_distr(j,i,:)))' * squeeze(c_prime(j,i,:));

    end
end

[mesh.k,mesh.m]=meshgrid(Q*grid.k,grid.m);

moneycapital_pdf=sum(joint_distr,3);

% Positive liquid only
% m_policy(grid.m(:)<0,:)=[];
% k_policy(grid.m(:)<0,:)=[];
% c_policy(grid.m(:)<0,:)=[];
% mesh.k(grid.m(:)<0,:)=[];
% moneycapital_pdf(grid.m(:)<0,:)=[];
% moneycapital_pdf=moneycapital_pdf./sum(moneycapital_pdf(:));
% mesh.m(grid.m(:)<0,:)=[];


switch(wealthsort)
    case('LiquidWealth')      
    mplusk=mesh.m(:);
    case('NetWealth')      
    mplusk=mesh.m(:)+mesh.k(:);
end
[mplusk, IX]=sort(mplusk);

mfull=mesh.m(IX);
kfull=mesh.k(IX);
Ifull=Inc_condh(IX)*par.W*par.N;

mplusIfull=mfull+1/6*Ifull;
m_policy=m_policy(IX);
k_policy=k_policy(IX);
c_policy=c_policy(IX);

moneycapital_pdf=moneycapital_pdf(IX);
moneycapital_cdf=cumsum(moneycapital_pdf);

%%
% Select percentiles
P=100;

% Select smoothing parameter
h=0.1;

lambda=ones(P,1);
lambdaLI=ones(P,1);
lambdaIL=ones(P,1);
lambdaC=ones(P,1);
lambdaLIP=ones(P,1);
lambdaILP=ones(P,1);
lambdaLIINC=ones(P,1);
locstdIL=ones(P,1);
locstdLI=ones(P,1);

data=[mplusk(:) mfull kfull moneycapital_cdf moneycapital_pdf*100 mplusIfull m_policy k_policy c_policy ];
data(data(:,5)==0,:)=[]; %drop combinations with zero prob
Y=(data(:,2))./(data(:,1)); % obtain li/il ratio
YLI=(data(:,2)); % Liquid assets
YIL=(data(:,3)); % Illiquid assets
YLIINC=(data(:,6)); % Adj liquid assets
YLIP=(data(:,7)); % Liquid assets policy
YILP=(data(:,8)); % Liquid assets policy
YC=(data(:,9)); % Consumption

for j=1:P
    
    
    dist=data(:,4)-j/P-0.5/P;
    weight=sqrt(normpdf(dist/h)).*sqrt(data(:,5));
    X=[weight dist.*weight];
   
    beta=regress(Y.*weight,X);
    lambda(j)=beta(1);
    beta=regress(YLI.*weight,X);
    lambdaLI(j)=(beta(1));
    beta=regress(YIL.*weight,X);
    lambdaIL(j)=(beta(1));
    beta=regress(YLIINC.*weight,X);
    lambdaLIINC(j)=(beta(1));
    beta=regress(YLIP.*weight,X);
    lambdaLIP(j)=(beta(1));
    beta=regress(YILP.*weight,X);
    lambdaILP(j)=(beta(1));
    beta=regress(YC.*weight,X);
    lambdaC(j)=(beta(1));

    weight=((abs(dist))<.2).*(data(:,5));
    LImean=sum(YLI.*weight)/sum(weight);
    ILmean=sum(YIL.*weight)/sum(weight);
    covvar=sum(((YLI-LImean)).*((YIL-ILmean).*weight))/sum(weight);
    locvarLI=sqrt((sum(((YLI-LImean).^2).*weight))/sum(weight))./LImean;
    locvarIL=sqrt((sum(((YIL-ILmean).^2).*weight))/sum(weight))./ILmean;
    locstdIL(j)=locvarIL;%covvar./sqrt(locvarIL*locvarLI);
    locstdLI(j)=locvarLI;
end


end
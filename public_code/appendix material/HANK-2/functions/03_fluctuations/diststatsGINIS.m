function [ginibonds,ginicapital,giniearnings,giniwealth,giniconsumption,giniincome] = diststatsGINIS(joint_distr,grid,mpar,c,R,RB,PI,W,naminus,nnminus,PROFITS,meshes,par)
% function [quints,gini,constrained] = diststats(joint_distr,kgrid)
% distats computes the following summary statistics given a
% distribution of agents over productivity and assets:
% - asset quintiles
% - gini coefficient
% - share of borrowing constrained agents
%
% author: Ralph L�tticke
% Date of last change: 04.11.2013

%% Money
gini.BC=sum(sum(joint_distr(1,:,:)));
% collapse distribution
money_pdf=sum(sum(joint_distr,2),3)';

% compute gini (using formular for discrete distributions, see
% http://en.wikipedia.org/wiki/Gini_coefficient)
S=cumsum(money_pdf.*grid.m);
S=[0 S];
ginibonds=1-(sum(money_pdf.*(S(1:end-1)+S(2:end)))/S(end));

%% Capital
capital_pdf=sum(sum(joint_distr,1),3);

% compute gini (using formular for discrete distributions, see
% http://en.wikipedia.org/wiki/Gini_coefficient)
S=cumsum(capital_pdf.*grid.k);
S=[0 S];
ginicapital=1-(sum(capital_pdf.*(S(1:end-1)+S(2:end)))/S(end));

%% Income
WW=((1-par.nu)*nnminus+par.nu*naminus).*W.*meshes.h/par.H;
WW(:,:,end)=PROFITS*par.profitshare;

[WW, IX]=sort(WW(:));

inc_pdf=joint_distr(IX);

S=cumsum(inc_pdf.*WW);
S=[0 S'];
giniearnings=1-(sum(inc_pdf'.*(S(1:end-1)+S(2:end)))/S(end));
% 
% h_pdf=squeeze(sum(sum(joint_distr,1),2))';
% S=cumsum(h_pdf.*W.*grid.h);
% S=[0 S];
% gini.earnings=1-(sum(h_pdf.*(S(1:end-1)+S(2:end)))/S(end));

%% M+K
for k = 1:mpar.nk
    for m = 1:mpar.nm
        
        mplusk(m+(k-1)*mpar.nm) = grid.m(m)+grid.k(k);
        
    end
end

[mplusk, IX]=sort(mplusk);

moneycapital_pdf=sum(joint_distr,3);
moneycapital_pdf=moneycapital_pdf(IX);
moneycapital_cdf=cumsum(moneycapital_pdf);
    
S=cumsum(moneycapital_pdf.*mplusk);
S=[0 S];
giniwealth=1-(sum(moneycapital_pdf.*(S(1:end-1)+S(2:end)))/S(end));

%% Total Income
WW=((1-par.nu)*nnminus+par.nu*naminus).*W.*meshes.h/par.H;
WW(:,:,end)=PROFITS*par.profitshare;

income=meshes.m*RB/PI+meshes.k*R+WW;

[income, IX]=sort(income(:));

inc_pdf=joint_distr(IX);

S=cumsum(inc_pdf.*income);
S=[0 S'];
giniincome=1-(sum(inc_pdf'.*(S(1:end-1)+S(2:end)))/S(end));

%% Consumption
% x=c;
% aux_x=par.tau*W*NH*meshes.h / (1+par.gamma);
% aux_x(:,:,end)=0;
% c=x+aux_x(:);

[c, IX]=sort(c(:));

c_pdf=joint_distr(IX);

S=cumsum(c_pdf.*c);
S=[0 S'];
giniconsumption=1-(sum(c_pdf'.*(S(1:end-1)+S(2:end)))/S(end));

end
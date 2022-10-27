% Steady state statistics
%%

clear targets
targets.ShareBorrower=sum((grid.m<0).*sum(joint_distr,2)');
targets.K=grid.m*(sum(joint_distr,2))-grid.B;
targets.B=grid.B;
targets.KY=targets.K/Output;
targets.BK=grid.B/targets.K;
targets.Y=Output;
BCaux_M=sum(sum(joint_distr,2),3);
targets.m_bc=BCaux_M(1,:);
targets.m_0=BCaux_M(grid.m==0);
par.RB = 1+R_fc;

labortax =(1-par.tau1).*W_fc.*N +(1-par.tau1).*Profits_fc;
par.gamma1=targets.B*(1-par.RB)+labortax(1);
par.W=W_fc(1);
par.PROFITS=Profits_fc(1);
par.N=N;
par.R=R_fc;
targets.GtoY=par.gamma1/Output;
targets.T=labortax;

targets.C = joint_distr(:)'*c_guess(:);%  B - Bminus*RBminus/PIminus + Tminus;

% Liquid Gini
[liquid_sort, IX]  = sort(meshes.m(:));
liquid_pdf         = joint_distr(IX);
liquid_cdf         = cumsum(liquid_pdf);
targets.Negliquid  = sum((liquid_sort<0).*liquid_pdf);

S                  = cumsum(liquid_pdf.*liquid_sort)';
S                  = [0 S];
targets.GiniLI      = 1-(sum(liquid_pdf.*(S(1:end-1)+S(2:end))')/S(end));

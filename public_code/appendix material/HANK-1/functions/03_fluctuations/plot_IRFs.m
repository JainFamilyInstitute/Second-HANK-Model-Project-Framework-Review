close all
mpar.maxlag=40;

IRF_distr=Gamma_state*IRF_state_sparse(1:numstates-4,1:mpar.maxlag);

IRF_distr=reshape(IRF_distr,[mpar.nm mpar.nh mpar.maxlag]);
IRF_distr_m=squeeze(sum(IRF_distr,2));
IRF_distr_h=squeeze(sum(IRF_distr,1));

% preparation
IRF_LIQUID_H=100*grid.h(1:end-1)*IRF_distr_h(1:end-1,2:end)/par.H;
K=grid.m*IRF_distr_m(:,1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_LIQUID_I=100*(I/(par.delta*grid.K)-1);
IRF_LIQUID_K=100*(K(2:end)./grid.K-1);

IRF_LIQUID_S=100*IRF_state_sparse(numstates-os+3,1:end-1);
BtoK=IRF_state_sparse(numstates-os+2,1:end)+grid.B/grid.K;

M=K.*BtoK;
IRF_LIQUID_M=100*(M(2:end)/grid.B-1);
IRF_LIQUID_Y=100*IRF_state_sparse(end-oc+3,1:end-1);
IRF_LIQUID_G=100*IRF_state_sparse(end-oc+4,1:end-1);
IRF_LIQUID_C=100*IRF_state_sparse(end-oc+12,1:end-1);
IRF_LIQUID_N=100*IRF_state_sparse(end-oc+8,1:end-1);
N=exp(IRF_state_sparse(end-oc+8,1:end))*par.N;
MC=exp(IRF_state_sparse(end-oc+13,1:end))*par.mu;

IRF_LIQUID_R=100*IRF_state_sparse(end-oc+6,1:end-1);
IRF_LIQUID_PI=100*100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_LIQUID_PROFITS=100*IRF_state_sparse(end-oc+7,1:end-1);
IRF_LIQUID_W=100*IRF_state_sparse(end-oc+5,1:end-1);

PI=1+IRF_state_sparse(end-oc+2,1:end-1);
Q=par.Q*(1+IRF_state_sparse(end-oc+1,1:end-1));
R=par.R*(1+IRF_state_sparse(end-oc+6,1:end-1));
Q=exp(log(par.Q)+IRF_state_sparse(end-oc+1,1:end-1));
R=exp(log(par.R)+IRF_state_sparse(end-oc+6,1:end-1));
RB=exp(log(par.RB)+(IRF_state_sparse(n2(2)+1,2:end)));

IRF_LIQUID_RB=100*100*(RB-par.RB);
IRF_LIQUID_RBREAL=100*100*(RB(1:end-1)./PI(2:end)-par.RB);
IRF_LIQUID_Q=100*100*(Q-par.Q);
IRF_LIQUID_D=100*100*((1+IRF_LIQUID_R/100)*par.R-par.R);
MK_Liquid=grid.B/grid.K;
M_Liquid=M;
K_Liquid=K;

IRF_RQ = 100*100*(((Q(2:end)+R(2:end))./Q(1:end-1)-(1+par.R)));

IRF_LIQUID_LP=100*100*(((Q(2:end)+R(2:end))./Q(1:end-1)-RB(1:end-1)./PI(2:end))-((1+par.R/par.Q)-par.RB));

save(['../saves/' casename '_IRF_LIQUID_data'],'IRF_LIQUID_Y','IRF_LIQUID_C','IRF_LIQUID_I','IRF_LIQUID_G','IRF_LIQUID_K','IRF_LIQUID_M','IRF_LIQUID_N','IRF_LIQUID_S','IRF_LIQUID_PI','IRF_LIQUID_RB','IRF_LIQUID_RBREAL','IRF_LIQUID_Q','IRF_LIQUID_D','IRF_LIQUID_LP','IRF_LIQUID_W','IRF_LIQUID_PROFITS','K_Liquid','M_Liquid','MK_Liquid') 


%%
close all



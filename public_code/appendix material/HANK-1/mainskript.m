
%% Initialize workspace and load directories
addpath(genpath('functions'))


%% Switch options
FindNewSS=true;
casename='SS_BASELINE';

mpar.overrideEigen = true;

%% Solve for Steady state
if FindNewSS
    
    % Set parameters
    defineSS_pars_BASELINE
    
    mainskript_steadystate
        
else
    load(casename)
    mpar.overrideEigen = true;

end

%% Select aggregate shock
aggrshock     = 'TFP';
par.sigmaS2=0.01;
par.rhoS2=0.95;


aggrshock     = 'MP';
par.sigmaS1=0.0009;
par.rhoS1=0.0;

[~,~,grid.boundsH] = Tauchen(par.rhoH,mpar.nh-1,1, 0, 'importance'); % LR variance = 1

grid.m=grid.m/(1+grid.B/grid.K);

casename='SS_BASELINE_PHI0';

par.phi=0;

%% Produce matrices to reduce state-space
mpar.maxdim=12;

mainskript_statereduc

disp('Computing system for SS.');
toc
F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets,Copula,P_H,aggrshock,oc);
tic
[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);
toc

%% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F, numstates,numcontrols,Gamma_state,oc,mpar,par,grid);

%% Produce IRFs to monetary shock
x0=zeros(numstates,1);
x0(end-1)=par.sigmaS1; % Monetary shock
% x0(end)=par.sigmaS2; % TFP shock

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;
mpar.maxlag=40;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

plot_IRFs
%%
casename='SS_BASELINE_PHI18';

par.phi=18;

%% Produce matrices to reduce state-space
mpar.maxdim=12;

mainskript_statereduc

disp('Computing system for SS.');
toc
F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets,Copula,P_H,aggrshock,oc);
tic
[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);
toc

%% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F, numstates,numcontrols,Gamma_state,oc,mpar,par,grid);

%% Produce IRFs to monetary shock
x0=zeros(numstates,1);
x0(end-1)=par.sigmaS1; % Monetary shock
% x0(end)=par.sigmaS2; % TFP shock

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;
mpar.maxlag=40;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

plot_IRFs


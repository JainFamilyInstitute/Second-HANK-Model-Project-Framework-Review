
%% Initialize workspace and load directories
addpath(genpath('functions'))

%% Select options
% Search for steady state and name it
FindNewSS           = true;
casename            = 'SS_BASELINE'


% Options
mpar.overrideEigen  = true; % Warning appears, but critical Eigenvalue shifted

%% Solve for Steady state
if FindNewSS
    % Set parameters
    defineSS_pars
    
    mainskript_steadystate
    
else
    load(casename)
    % Options
    mpar.overrideEigen  = true; % Warning appears, but critical Eigenvalue shifted
end

%% Select aggregate shock
aggrshock           = 'MP';
par.rhoS            = 0.0;    % Persistence of variance
par.sigmaS          = 0.0009 ;    % quarter STD of variance shocks

%% Produce matrices to reduce state-space
mpar.maxdim=12;

mainskript_statereduc

%% Initialize System

disp('Computing system for SS.');
F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets,Copula,P_H,aggrshock,oc);
tic
[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);
toc

%% Solve RE via Schmitt-Grohe-Uribe Form
fprintf('Taylor Rule Case: %d\n',TayRule_case);
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F, numstates,numcontrols,Gamma_state,oc,mpar,par,grid);

%% Produce IRFs (Figures 1, 2, 3)
mpar.maxlag= 40; % Quarters

x0=zeros(numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

casename = save_IRF_data;
plot_IRFs_prep


%% Produce Table 3
redistribution_table

%% General Equilibrium Liquidity (Figure 7)
mpar.maxlag= 40; % Quarters

x0=zeros(numstates,1);
x0(end)=0.0024; % average annual shock

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

quart=12;
[~,~,~,JD_ss,c_a_star_ss,m_a_star_ss,k_a_star_ss,c_n_star_ss,m_n_star_ss]=F(State,State_m,Contr,Contr_m);
[~,~,~,JD_a1,c_a_star_a1,m_a_star_a1,k_a_star_a1,c_n_star_a1,m_n_star_a1]=F(IRF_state_sparse(1:numstates,quart-1),IRF_state_sparse(1:numstates,quart-2),IRF_state_sparse(numstates+1:end,quart-1),IRF_state_sparse(numstates+1:end,quart-2));
[~,~,~,JD_a2,c_a_star_a2,m_a_star_a2,k_a_star_a2,c_n_star_a2,m_n_star_a2]=F(IRF_state_sparse(1:numstates,quart),IRF_state_sparse(1:numstates,quart-1),IRF_state_sparse(numstates+1:end,quart),IRF_state_sparse(numstates+1:end,quart-1));

portfolio_response



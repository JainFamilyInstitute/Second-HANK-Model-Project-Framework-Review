function [Difference,LHS,RHS,JD_new,c_a_star,m_a_star,k_a_star,c_n_star,m_n_star,P]  = Fsys(State,Stateminus,...
    Control_sparse,Controlminus_sparse,StateSS,...
    ControlSS,Gamma_state,Gamma_control,InvGamma,...
    par,mpar,grid,targets,Copula,P,aggrshock,oc)
% System of equations written in Schmitt-Grohe-Uribe generic form with states and controls
% STATE: Vector of state variables t+1 (only marginal distributions for histogram)
% STATEMINUS: Vector of state variables t (only marginal distributions for histogram)
% CONTROL: Vector of state variables t+1 (only coefficients of sparse polynomial)
% CONTROLMINUS: Vector of state variables t (only coefficients of sparse polynomial)
% STATESS and CONTROLSS: Value of the state and control variables in steady
% state. For the Value functions these are at full grids.
% GAMMA_STATE: Mapping such that perturbationof marginals are still
% distributions (sum to 1).
% GAMMA_CONTROL: Values of the polynomial base at all nodes to map
% sparse coefficient changes to full grid
% INVGAMMA: Projection of Value functions etc. to Coeffeicent space for
% sparse polynomials.
% PAR, MPAR: Model and numerical parameters (structure)
% GRID: Liquid, illiquid and productivity grid
% TARGETS: Stores targets for government policy
% COPULA: Interpolant that allows to map marginals back to full-grid
% distribuitions
% P: steady state transition matrix
% aggrshock: sets wether the Aggregate shock is TFP or uncertainty
%
%% Initializations
mutil = @(c)(1./(c.^par.xi));
invmutil = @(mu)((1./mu).^(1/par.xi));

% Number of states, controls
nx   = mpar.nm+mpar.nk+mpar.nh -4 +2; % Number of states
ny   = length(ControlSS); % number of Controls
NxNx = nx-2; % Number of states without aggregate
NN   = mpar.nm*mpar.nh*mpar.nk; % Number of points in the full grid

% Initialize LHS and RHS
LHS  = zeros(nx+ny,1);
RHS  = zeros(nx+ny,1);

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_a_cind = (1:NN);
mutil_n_cind = 1*NN+(1:NN);
Vkind = 2*NN+(1:NN);
Qind  = 3*NN+1;
PIind = 3*NN+2;
Yind  = 3*NN+3;
Gind  = 3*NN+4;
Wind  = 3*NN+5;
Rind  = 3*NN+6;
Profitind  = 3*NN+7;
Nind  = 3*NN+8;
Tind  = 3*NN+9;
Kind  = 3*NN+10;
Bind  = 3*NN+11;
Cind  = 3*NN+12;
MCind  = 3*NN+13;

% Indexes for States
marginal_mind = (1:mpar.nm-1);
marginal_kind = (mpar.nm-1+(1:mpar.nk-1));
marginal_hind = (mpar.nm+mpar.nk-2 + (1:(mpar.nh-2)));

RBind = NxNx+1;
Sind  = NxNx+2;

%% Control Variables (Change Value functions according to sparse polynomial)
Control      = ControlSS .* (1+Gamma_control*(Control_sparse));
Controlminus = ControlSS .* (1+Gamma_control*(Controlminus_sparse));

Control(end-oc+1:end)       = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Control_sparse);
Controlminus(end-oc+1:end)  = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Controlminus_sparse);


%% State Variables
% read out marginal histogramm in t+1, t
Distribution      = StateSS(1:end-2) + Gamma_state * State(1:NxNx);
Distributionminus = StateSS(1:end-2) + Gamma_state * Stateminus(1:NxNx);

% Aggregate Endogenous States
RB      = StateSS(end-1) + (State(end-1));
RBminus = StateSS(end-1) + (Stateminus(end-1));

% Aggregate Exogenous States
S       = StateSS(end) + (State(end));
Sminus  = StateSS(end) + (Stateminus(end));

%% Split the Control vector into items with names
% Controls
Vk       = mutil(Control(Vkind));
Vkminus  = mutil(Controlminus(Vkind));
mutil_c_a       = mutil(Control(mutil_a_cind));
mutil_c_a_minus  = mutil(Controlminus(mutil_a_cind));
mutil_c_n       = mutil(Control(mutil_n_cind));
mutil_c_n_minus  = mutil(Controlminus(mutil_n_cind));

% Aggregate Controls (t+1)
PI = exp(Control(PIind));
Y  = exp(Control(Yind ));
K  = exp(Control(Kind ));
B  = exp(Control(Bind ));

% Aggregate Controls (t)
PIminus = exp(Controlminus(PIind));
Qminus  = exp(Controlminus(Qind ));
Yminus  = exp(Controlminus(Yind ));
Gminus  = exp(Controlminus(Gind ));
Wminus  = exp(Controlminus(Wind ));
Rminus  = exp(Controlminus(Rind ));
Profitminus  = exp(Controlminus(Profitind ));
Nminus  = exp(Controlminus(Nind ));
Tminus  = exp(Controlminus(Tind ));
Kminus  = exp(Controlminus(Kind ));
Bminus  = exp(Controlminus(Bind ));
Xminus  = exp(Controlminus(Cind ));
mcminus  = exp(Controlminus(MCind ));

%% Write LHS values
% Controls
LHS(nx+Vkind)      = invmutil(Vkminus);
LHS(nx+mutil_a_cind) = invmutil(mutil_c_a_minus);
LHS(nx+mutil_n_cind) = invmutil(mutil_c_n_minus);

LHS(nx+Qind)       = (Qminus);
LHS(nx+Yind)       = (Yminus);

LHS(nx+Gind)       = (Gminus);
LHS(nx+Wind)       = (Wminus);
LHS(nx+Rind)       = (Rminus);
LHS(nx+Profitind)  = (Profitminus);
LHS(nx+Nind)       = (Nminus);
LHS(nx+Tind)       = (Tminus);
LHS(nx+Kind)       = (Kminus);
LHS(nx+Bind)       = (Bminus);
LHS(nx+Cind)       = (Xminus);
LHS(nx+MCind)      = (mcminus);

% States
% Marginal Distributions (Marginal histograms)
LHS(marginal_mind) = Distribution(1:mpar.nm-1);
ba=mpar.nm;
LHS(marginal_kind) = Distribution(ba+(1:mpar.nk-1));
ba=mpar.nm+mpar.nk;
LHS(marginal_hind) = Distribution(ba+(1:mpar.nh-2));
LHS(RBind)         = (RB);
LHS(Sind)          = (S);

% take into account that RB is in logs
RB=exp(RB); RBminus=exp(RBminus);
Tauminus = par.tau;

%% Set of Differences for exogenous process
RHS(Sind) = (par.rhoS * (Sminus));

switch(aggrshock)
    case('MP')
        EPS_TAYLOR=-Sminus;
        TFP=1;
        
    case('TFP')
        TFP=exp(Sminus);
        EPS_TAYLOR=0;
end

% Calculate aggregate Capital, Bonds and Human Capital Supply in t
marginal_mminus = Distributionminus(1:mpar.nm)';
marginal_kminus = Distributionminus((1:mpar.nk)+mpar.nm)';
marginal_hminus = Distributionminus(mpar.nm+mpar.nk+(1:mpar.nh))';

Hminus  = sum(grid.h(1:end-1).*marginal_hminus(1:end-1)); %Last column is entrepreneurs.
Lminus  = sum(grid.m.*marginal_mminus);

RHS(nx+Bind) = Lminus ;
RHS(nx+Kind) = (sum(grid.k.*marginal_kminus));

% Calculate joint distributions
cumdist = zeros(mpar.nm+1,mpar.nk+1,mpar.nh+1);
cumdist(2:end,2:end,2:end) = Copula({cumsum(marginal_mminus),cumsum(marginal_kminus)',cumsum(marginal_hminus)});
JDminus = diff(diff(diff(cumdist,1,1),1,2),1,3);

% Generate meshes for b,k,h
[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

%% Aggregate Output
mc              =  par.mu- (par.beta * log(PI)*Y/Yminus - log(PIminus))/par.kappa;
mutil_c         =  par.nu*mutil_c_a(:)+(1-par.nu)*mutil_c_n;

% Incomes (grids)
inc.rent    = meshes.k*Rminus;
inc.capital = meshes.k*Qminus;
inc.money   = (RBminus/PIminus).*meshes.m...
    + (meshes.m<0).*(par.borrwedge/PIminus).*meshes.m;


%% First Set: Value Functions, Marginal Values
%% Update policies

EVk = reshape(reshape(Vk,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);
RBaux = RB/PI + (meshes.m<0).*(par.borrwedge/PI);
EVm = reshape(reshape(RBaux(:).*mutil_c,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);

[c_a_star,m_a_star,k_a_star,c_n_star,m_n_star] = EGM_policyupdate(EVm,EVk,Qminus,PIminus,RBminus,Wminus,Profitminus,Tauminus,TFP*mc,Kminus,Nminus,inc,meshes,grid,par,mpar);

naminus= ((1/par.gamma).* mutil(c_a_star(:)).*meshes.h(:).*par.tau.*par.alpha.*TFP.*mc.*(Kminus./Nminus).^(1-par.alpha)).^(par.sigma_n);
naminus          = min(naminus,1);
naminus          = reshape(naminus,[mpar.nm mpar.nk mpar.nh]);
naminus(:,:,end) = 0;
nnminus= ((1/par.gamma).* mutil(c_n_star(:)).*meshes.h(:).*par.tau.*par.alpha.*TFP.*mc.*(Kminus./Nminus).^(1-par.alpha)).^(par.sigma_n);
nnminus          = min(nnminus,1);
nnminus          = reshape(nnminus,[mpar.nm mpar.nk mpar.nh]);
nnminus(:,:,end) = 0;

N  = (par.nu.*JDminus(:)'*naminus(:)+(1-par.nu).*JDminus(:)'*nnminus(:));


%% Prices that are not part of Control Vector
RHS(nx+Nind)    = N;
RHS(nx+Yind)    = (TFP*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));

% Wage Rate
RHS(nx+Wind) =  TFP *par.alpha       *mc.* (Kminus./(Nminus)).^(1-par.alpha);
% Return on Capital
RHS(nx+Rind) =  TFP *(1-par.alpha)   *mc.* ((Nminus)./Kminus).^(par.alpha)  - par.delta;

% Profits for Entrepreneurs
RHS(nx+Profitind) = (1-mc)*Yminus - Yminus.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2 + 1/2*par.phi*((K-Kminus).^2)./Kminus;

meshaux=meshes;
meshaux.h(:,:,end)=1000;

%% Update Marginal Value Bonds
mutil_c_n = mutil(c_n_star); % marginal utility at consumption policy no adjustment
mutil_c_a = mutil(c_a_star); % marginal utility at consumption policy adjustment
% mutil_c_aux    = par.nu.*mutil_c_a + (1-par.nu).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)
RHS(nx+mutil_a_cind) = invmutil(mutil_c_a(:)); % Write Marginal Utility to RHS of F
RHS(nx+mutil_n_cind) = invmutil(mutil_c_n(:)); % Write Marginal Utility to RHS of F
%% Update marginal Value of Capital
EVk     = reshape(Vk,[mpar.nm*mpar.nk, mpar.nh])*P';
Vk_next = griddedInterpolant(meshaux.m,meshaux.k,meshaux.h,reshape(EVk,[mpar.nm mpar.nk mpar.nh]));
Vk_aux  = par.nu.*(Rminus+Qminus).*mutil_c_a + (1-par.nu).*Rminus.*mutil_c_n+ par.beta.*(1-par.nu).*Vk_next(m_n_star,meshaux.k,meshaux.h); % Expected marginal utility at consumption policy (w &w/o adjustment)
RHS(nx+Vkind) = invmutil(Vk_aux(:));  % Write Marginal Value of K to RHS of F


%% Differences for distributions

% Transitions: find next smallest on-grid value for money and capital choices
[Dist_m_a,idm_a] = genweight(m_a_star,grid.m);
[Dist_m_n,idm_n] = genweight(m_n_star,grid.m);
[Dist_k,idk_a]   = genweight(k_a_star,grid.k);

% Iterate Distribution function forward
JD_new = JDIteration_mex(reshape(JDminus,[mpar.nm, mpar.nk, mpar.nh]),...
    int32(idm_n),int32(idm_a),int32(idk_a),Dist_m_n,Dist_m_a,Dist_k,...
    par.nu*ones([mpar.nm, mpar.nk, mpar.nh]),P);

JD_new = reshape(JD_new(:),[mpar.nm,mpar.nk,mpar.nh]);

% Next period marginal histograms
% liquid assets
aux_m = squeeze(sum(sum(JD_new,2),3));
RHS(marginal_mind) = aux_m(1:end-1); %Leave out last state
% illiquid assets
aux_k =  squeeze(sum(sum(JD_new,1),3))' ;
RHS(marginal_kind) = aux_k(1:end-1); %Leave out last state
% human capital
aux_h = squeeze(sum(sum(JD_new,1),2));
RHS(marginal_hind) = aux_h(1:end-2); %Leave out last state & entrepreneurs

%% Third Set: Government Budget constraint
% Return on bonds (Taylor Rule)
% RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR;
global TayRule_case

if TayRule_case == 1
    RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR-0.01;

elseif TayRule_case == 2
    RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR...
       -0.01 + log(N/targets.L)*(1-par.rho_R); 

else
    RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR...
       -0.01 + min(log(N/targets.L),0)*(1-par.rho_R)*0.85; 
end

TN=(1-Tauminus).*(nnminus).*Wminus.*meshes.h./par.H;
TA=(1-Tauminus).*(naminus).*Wminus.*meshes.h./par.H;
TT=[TA(:); TN(:)];
TT=[par.nu.*JDminus(:);(1-par.nu).*JDminus(:)]'*TT(:);
taxrevenue =TT +(1-Tauminus).*Profitminus;

% Inflation jumps to equilibrate real bond supply and demand
RHS(nx+PIind) = par.rho_B * log((Bminus)/(targets.B)) ...
    + par.rho_B * log(RBminus/par.RB)...
    - (par.rho_B+par.gamma_pi) * log(PIminus/par.PI) ...
    - par.gamma_T * log((Tminus)/(targets.T));
LHS(nx+PIind) = log((B)/(targets.B));

% Government expenditures
RHS(nx+Gind) = B - Bminus*RBminus/PIminus + Tminus;

RHS(nx+Tind) = taxrevenue;

RHS(nx+MCind) = 1-mc;


% Resulting Price of Capital
RHS(nx+Qind)=(par.phi*(K./Kminus-1)+1)-par.ABS;

x=par.nu.*JDminus(:)'*c_a_star(:)+(1-par.nu).*JDminus(:)'*c_n_star(:);
RHS(nx+Cind) = x;
%% Difference
Difference=InvGamma(:,:)*((LHS-RHS)./[ones(nx,1);(ControlSS(1:end-oc));ones(oc,1)]);

end
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function

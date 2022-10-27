function [Difference,LHS,RHS,JD_new,c_star,m_star,P]  = Fsys(State,Stateminus,...
    Control_sparse,Controlminus_sparse,StateSS,...
    ControlSS,Gamma_state,Gamma_control,InvGamma,...
    par,mpar,grid,targets,Copula,P,aggrshock,oc)
% System of equations written in Schmitt-Grohé-Uribe generic form with states and controls
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


%% Initializations
util  = @(c)(c.^(1-par.sigma))./(1-par.sigma);
mutil = @(c)(1./(c.^par.sigma));
invutil = @(u)(((1-par.sigma).*u).^(1/(1-par.sigma)));
invmutil = @(mu)((1./mu).^(1/par.sigma));
% Number of states, controls
nx   = mpar.nm*mpar.nh +4 - mpar.nh -1; % Number of states
ny   = length(ControlSS); % number of Controls
NxNx = nx-4; % Number of states without aggregate
NN   = mpar.nm*mpar.nh; % Number of points in the full grid

% Initialize LHS and RHS
LHS  = zeros(nx+ny,1);
RHS  = zeros(nx+ny,1);

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_cind = 1:NN;
Qind  = 1*NN+1;
PIind = 1*NN+2;
Yind  = 1*NN+3;
Gind  = 1*NN+4;
Wind  = 1*NN+5;
Rind  = 1*NN+6;
Profitind  = 1*NN+7;
Nind  = 1*NN+8;
Tind  = 1*NN+9;
Kind  = 1*NN+10;
Bind  = 1*NN+11;
Cind  = 1*NN+12;
MUind  = 1*NN+13;

% Indexes for States
distr_ind = 1:mpar.nm*mpar.nh-mpar.nh-1;

RBind = NxNx+1;
BtoKind  = NxNx+2;
S1ind  = NxNx+3;
S2ind  = NxNx+4;

%% Control Variables (Change Value functions according to sparse polynomial)
Control      = ControlSS+Control_sparse;
Controlminus = ControlSS+Controlminus_sparse;

%% State Variables
% read out marginal histogramm in t+1, t
Distribution      = StateSS(1:end-4) + Gamma_state * State(1:NxNx);
Distributionminus = StateSS(1:end-4) + Gamma_state * Stateminus(1:NxNx);

% Aggregate Endogenous States
RB      = StateSS(end-3) + (State(end-3));
RBminus = StateSS(end-3) + (Stateminus(end-3));

% Aggregate Exogenous States
BtoK       = StateSS(end-2) + (State(end-2));
BtoKminus  = StateSS(end-2) + (Stateminus(end-2));

% Aggregate Exogenous States
S1       = StateSS(end-1) + (State(end-1));
S1minus  = StateSS(end-1) + (Stateminus(end-1));
S2       = StateSS(end) + (State(end));
S2minus  = StateSS(end) + (Stateminus(end));
%% Split the Control vector into items with names
% Controls
mutil_c       = mutil(Control(mutil_cind));
mutil_cminus  = mutil(Controlminus(mutil_cind));

% Aggregate Controls (t+1)
PI = exp(Control(PIind));
Y  = exp(Control(Yind ));
K  = exp(Control(Kind ));
Q  = exp(Control(Qind ));
R  = exp(Control(Rind ));
B = exp(Control(Bind ));

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
Cminus  = exp(Controlminus(Cind ));
MUminus  = exp(Controlminus(MUind ));

%% Write LHS values
% Controls
LHS(nx+mutil_cind) = invmutil(mutil_cminus);

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
LHS(nx+Cind)     = (Cminus);
LHS(nx+MUind)      = (MUminus);

% States
% Marginal Distributions (Marginal histograms)
LHS(distr_ind) = Distribution(1:mpar.nm*mpar.nh-1-mpar.nh);

LHS(RBind)         = (RB);
LHS(S1ind)          = (S1);
LHS(S2ind)          = (S2);

LHS(BtoKind)       = (BtoK);

Tauminus=par.tau1;

% take into account that RB and B are in logs
RB=exp(RB); RBminus=exp(RBminus);
BtoK=exp(BtoK); BtoKminus=exp(BtoKminus);

%PI=exp(PI); PIminus=exp(PIminus);
%% Set of Differences for exogenous process
RHS(S1ind) = (par.rhoS1 * (S1minus));
RHS(S2ind) = (par.rhoS2 * (S2minus));

        TFP=exp(S2minus);
        EPS_TAYLOR=-S1minus;

marginal_mminus = sum(reshape(Distributionminus,[mpar.nm mpar.nh]),2);
marginal_hminus = squeeze(sum(reshape(Distributionminus,[mpar.nm mpar.nh]),1));

Hminus  = grid.h(1:end-1)*marginal_hminus(1:end-1)'; %Last column is entrepreneurs.

RHS(nx+Kind)= grid.m*marginal_mminus;

% Generate meshes for b,k,h
[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

RHS(nx+MUind)             =  par.mu- (par.beta * log(PI)*Y/Yminus - log(PIminus))/par.kappa;


% Incomes (grids)
inc.money   = (Rminus+Qminus+BtoKminus*RBminus/PIminus)*meshes.m;

%% Update policies

Raux = (R+Q+BtoK*RB/PI)/(Qminus+BtoK); % Expected marginal utility at consumption policy (w &w/o adjustment)
EVm = reshape(reshape(Raux(:).*mutil_c,[mpar.nm mpar.nh])*P',[mpar.nm, mpar.nh]);

[c_star,m_star] = EGM_policyupdate(EVm,Rminus,Qminus,BtoK,BtoKminus,RBminus,PIminus,Hminus,Wminus,Kminus,Nminus,TFP*MUminus,Profitminus,Tauminus,inc,meshes,grid,par,mpar);

n_star= ((1/par.gamma).*mutil(c_star(:)).*meshes.h(:).*Tauminus.*Wminus).^(par.sigma_n);
n_star          = min(n_star,1);
n_star          = reshape(n_star,[mpar.nm mpar.nh]);
n_star(:,end)   = 0;

Nguess = Distributionminus(:)'*n_star(:);    

RHS(nx+Nind)    = Nguess;
RHS(nx+Yind)    = (TFP*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));

% Prices that are not part of Control Vector

% Wage Rate
RHS(nx+Wind) = TFP *par.alpha       *MUminus.* (Kminus./(Nminus)).^(1-par.alpha);
% Return on Capital
RHS(nx+Rind) = TFP *(1-par.alpha)   *MUminus.* ((Nminus)./Kminus).^(par.alpha)  - par.delta;

% Profits for Entrepreneurs
RHS(nx+Profitind) = (1-MUminus)*Yminus - Yminus.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2 + 1/2*par.phi*((K-Kminus).^2)./Kminus;

%% Update Marginal Value Bonds
RHS(nx+mutil_cind) = c_star(:);%invmutil(mutil_c_aux(:)); % Write Marginal Utility to RHS of F

%% Differences for distributions

% find next smallest on-grid value for money and capital choices
weight11  = zeros(mpar.nm, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm, mpar.nh,mpar.nh);

% Adjustment case
[Dist_m,idm] = genweight(m_star,grid.m);

idm=repmat(idm(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nh],idm(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nh],idm(:)+1,idh(:));

for hh=1:mpar.nh
    
    %Corresponding weights
    weight11_aux = (1-Dist_m(:,hh));
    weight12_aux =  (Dist_m(:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P(hh,:);
end

weight11=permute(weight11,[1 3 2]);
weight12=permute(weight12,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nh,[1 2*mpar.nh]);

H=sparse(rowindex,[index11(:); index12(:)],...
    [weight11(:); weight12(:)],mpar.nm*mpar.nh,mpar.nm*mpar.nh); % mu'(h',k'), a without interest


JD_new=Distributionminus(:)'*H;

JD_new = reshape(JD_new(:),[mpar.nm,mpar.nh]);

% Next period marginal histograms
RHS(distr_ind) = JD_new(1:mpar.nm*mpar.nh-1-mpar.nh); %Leave out last state


%% Third Set: Government Budget constraint
% Return on bonds (Taylor Rule)
global TayRule_case
switch TayRule_case 

    case 1
    RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR;

    case 2
    RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR...
        + log(Nminus/targets.L).*(1-par.rho_R); % + log(Yminus/targets.Y).*(1-par.rho_R)*0.5 ;

    case 3
    RHS(RBind) = log(par.RB) + par.rho_R*log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi)+EPS_TAYLOR...
        + min(log(Nminus/targets.L),0).*(1-par.rho_R) ;
end

labortax =(1-Tauminus).*Wminus.*Nminus +(1-Tauminus).*Profitminus;

% Inflation jumps to equilibrate real bond supply and demand
RHS(nx+PIind) = par.rhoB * log((Bminus)/(targets.B)) ...
              + par.rhoB * log(RBminus/par.RB)...
              - (par.rhoB+par.gamma3) * log(PIminus) ...
              - par.gamma4 * log((Yminus)/(targets.Y));
LHS(nx+PIind) = log((B)/(targets.B));

RHS(nx+Bind) = BtoKminus*Kminus;

LHS(BtoKind) = RB/PI;
RHS(BtoKind) = R/Qminus+Q/Qminus;% log(B/K);


% Government expenditures
% RHS(nx+Gind) =  Gov_spend;
RHS(nx+Gind) = B - Bminus*RBminus/PIminus + Tminus;

RHS(nx+Tind) = labortax;

RHS(nx+Qind)=(par.phi*(K./Kminus-1)+1);

x=c_star;
RHS(nx+Cind) = Distributionminus(:)'*x(:);

%% Difference
Difference=(LHS-RHS);%InvGamma(:,:)*((LHS-RHS)./[ones(nx,1);(ControlSS(1:end-oc));ones(oc,1)]);


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

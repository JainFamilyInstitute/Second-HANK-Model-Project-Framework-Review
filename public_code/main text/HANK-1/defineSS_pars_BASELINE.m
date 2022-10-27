% Parameters
% Household Parameters
par.beta      = 0.984;      % discount factor
par.sigma     = 4;          % CRRA
par.gamma     = 1;
par.sigma_n   = 1;
par.probit2   = 0;
par.probitratio = 0;
par.probit1     = par.probitratio*par.probit2;
par.nu          = 1;
par.guessadj    = par.nu*1/(1+exp(par.probit1/par.probit2));

% Income Process
par.rhoH        = 0.98;    % persistence of productivity
par.sigmaH      = 0.06;    % STD of productivity shocks
par.rhoS        = 0.0;    % persistence of variance
par.sigmaS      = 0.0;    % STD of variance shocks
mpar.in         = 0.0006;  % prob. to become entrepreneur
mpar.out        = 0.0625;   % prob. to become worker

% Firm Side Parameters
par.mu          = 19/20;        % Markup 5%
par.alpha       = 2/3/par.mu;
par.delta       = 0.054/4;          % Depreciation rate
par.phi         = 10;

% Phillips Curve
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)

% Central Bank Policy
par.theta_pi    = 1.5;
par.rho_RB      = 0.8;

% Tax Schedule
par.tau0        = 0;         % lump-sum
par.tau1        = 0.7;         % proportional (net income)
par.tau2        = 0.0;         % progessivity
par.tauK        = 0.0;

%% Debt rule

par.rhoB=0.86;
par.gamma2 = 0;
par.gamma3 = 0;  
par.gamma4 = 0;

%% Returns
par.PI  = 1.02^.25; % inflation
par.EPI = 0; % expected inflation

par.RB  = (par.PI*1.02)^0.25; % 0; % real return times inflation

par.ABS = 0.0; % Loan to value ratio max.
par.mortwedge = 1.00^0.25-1; % wedge on mortgages 
par.borrwedge   = 1.16^0.25-1; % wedge on borrowing beyond secured borrowing
par.Q  = 1*(1-par.ABS);

%% Grids
% Idiosyncratic States
mpar.nm         = 150;
mpar.nk         = 0;
mpar.nh         = 4;

% Aggr States
mpar.nK         = 1;
mpar.nM         = 1;
mpar.ns         = 1;

%% Numerical Parameters
mpar.crit    = 1e-8;

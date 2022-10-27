function [c_a_guess,c_n_guess, psi_guess, inc]=policyguess(meshes,WWa,WWn,RR,RBRB,par,mpar)
%policyguess returns autarky policy guesses (in the first period only):
% c_a_guess, c_n_guess, psi_guess
% as well as income matrices later on used in EGM:
% inc
%
% Consumption is compositite leisure and physical consumption (x_it) in the
% paper, therefore labor income is reduced by the fraction of leisure
% consumed.

inc.labor_a   = par.tau.*WWa.*meshes.h;
inc.labor_n   = par.tau.*WWn.*meshes.h;

inc.rent    = RR.*meshes.k;
inc.money   = (par.RB/par.PI).*meshes.m...
    + (meshes.m<0).*(par.borrwedge/par.PI).*meshes.m;

inc.capital = par.Q.*meshes.k;
inc.profits = zeros(mpar.nm,mpar.nk,mpar.nh*mpar.ns);%lump sum profits

%% Initial policy guesses:  Autarky policies as guess

% Consumption guess
c_a_guess = inc.labor_a + inc.rent + inc.capital + max(inc.money,0) + inc.profits;
c_n_guess = inc.labor_n + inc.rent               + max(inc.money,0) + inc.profits;

% Initially guessed marginal continuation value of holding capital
psi_guess      = zeros([mpar.nm mpar.nk mpar.nh*mpar.ns]);

end

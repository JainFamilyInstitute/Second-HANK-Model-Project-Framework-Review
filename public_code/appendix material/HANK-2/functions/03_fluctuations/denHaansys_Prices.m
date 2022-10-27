function [JD,RB,EqQuant,PF,W,naminus,nnminus,R,c_a,m_a,k_a,c_n,m_n,P]  = denHaansys_Prices(JDminus,RBminus,Sminus,...
                            Control_F,Control_GUESS,ControlSS,Gamma_control,par,mpar,grid,P,aggrshock,oc,targets)
% denHaansys iterates the states forward using expected prices and value funtions
% that stored in the forecast CONTROL_F. It uses prices as stored in Control_GUESS
% as starting guesses and finds the within period market clearing prices.
%
% GAMMA_CONTROL maps polynomial coefficients to the value function (and its marginals),
% see Fsys. Parameters are stored in  par and mpar.
% JDminus, RBminus and Sminus contain the time t state variables.
% The function returns time t+1 state variables JD (distribution), RB (CB interest rate)
% and S (uncertainty state), along policy functions (c_a,m_a,k_a,c_n,m_n).
%
%
% Authors: Christian Bayer and Ralph C. Luetticke, Bonn and London, March
% 2017%
% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Precautionary Savings, Illiquid Assets, and the Aggregate Consequences of
%  Shocks to Household Income Risk', Bonn mimeo
% http://wiwi.uni-bonn.de/hump/wp.html
% =========================================================================
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The most recent version or successor of the above-mentioned paper is
% properly cited in all work and publications that benefit from insights
% derived from the Software.
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% =========================================================================


%% Initializations
util  = @(c)(c.^(1-par.xi))./(1-par.xi);
mutil = @(c)(1./(c.^par.xi));

% Number of states, controls
NN   = mpar.nm*mpar.nh*mpar.nk; % Number of points in the full grid


%% Indexes for LHS/RHS
% Indexes for Controls
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


%% Control Variables (Change Value functions according to sparse polynomial)
Control      = ControlSS .* (1+Gamma_control*(Control_F));
Controlminus = ControlSS .* (1+Gamma_control*(Control_GUESS));

Control(end-oc+1:end)       = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Control_F);
Controlminus(end-oc+1:end)  = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Control_GUESS);



%% Split the Control vector into items with names
% Controls
Vk       = mutil(Control(Vkind));
mutil_c_a       = mutil(Control(mutil_a_cind));
mutil_c_n       = mutil(Control(mutil_n_cind));


% Aggregate Controls (t+1)
PI = exp(Control(PIind));
Y  = exp(Control(Yind ));
% Aggregate Controls (t)
PIminus = exp(Controlminus(PIind));
Qminus  = exp(Controlminus(Qind ));
Yminus  = exp(Controlminus(Yind ));
Gminus  = exp(Controlminus(Gind ));
Wminus  = exp(Controlminus(Wind ));
Rminus  = exp(Controlminus(Rind ));
Profitminus  = exp(Controlminus(Profitind ));
Nminus  = exp(Controlminus(Nind ));


% States
% Marginal Distributions (Marginal histograms)

% Calculate aggregate Capital, Bonds and Human Capital Supply in t
marginal_hminus=squeeze(sum(sum(JDminus,1),2));
marginal_kminus=squeeze(sum(sum(JDminus,1),3));
marginal_bminus=squeeze(sum(sum(JDminus,2),3));

Hminus  = grid.h(1:end-1)*marginal_hminus(1:end-1); %Last column is entrepreneurs.
Kminus  = sum(grid.k.*marginal_kminus);
Bminus  = sum(grid.m.*marginal_bminus');

% take into account that RB is in logs
RBminus=exp(RBminus);

switch(aggrshock)
    case('Uncertainty')
        [P,~,~] = ExTransitions(exp(Sminus),grid,mpar,par);
end

DD = @(Qg)PriceD(Qg,Vk,PI,Y/Yminus,Bminus,Kminus,Nminus,Hminus,RBminus,Sminus,mutil_c_a,mutil_c_n,mpar,par,grid,aggrshock,P,JDminus,targets);
QuantGuess=[Qminus,PIminus];

options=optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

[EqQuant]=fsolve(@(p)DD(p),QuantGuess,options);
[~,JD,RB,PF,W,naminus,nnminus,R,c_a,m_a,k_a,c_n,m_n,AProb]= DD(EqQuant);
fprintf('%2.2f %8.2f \n',(abs(EqQuant-QuantGuess)./QuantGuess)*100);

end

function [QuantDiff,JD,RB,Profitminus,Wminus,naminus,nnminus,Rminus,c_a_star,m_a_star,k_a_star,c_n_star,m_n_star,AProb] = PriceD(QuantGuess,Vk,PI,YG,Bminus,Kminus,Nminus,Hminus,RBminus,Sminus,mutil_c_a,mutil_c_n,mpar,par,grid,aggrshock,P,JDminus,targets)
Qminus=QuantGuess(1);
PIminus=QuantGuess(2);

switch(aggrshock)
    case('TFP')
        TFP=exp(Sminus);
        MPshock=0;
    case('Uncertainty')
         % Tauchen style for Probability distribution next period
        [P,~,~] = ExTransitions(exp(Sminus),grid,mpar,par);
        TFP=1;%(Hminus/par.H).^(-par.alpha);
        MPshock=0;
    case('MP')
        TFP=1;
        MPshock=Sminus;
end

K  = ((Qminus-1+par.ABS)/par.phi+1)*Kminus;

RB = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi) + MPshock;

util  = @(c)(c.^(1-par.xi))./(1-par.xi);
mutil = @(c)(1./(c.^par.xi));

Tauminus  = par.tau;%*exp(-par.gamma_T*log(Bminus./(targets.B)));


mc        =  par.mu- (par.beta * log(PI)*YG - log(PIminus))/par.kappa;

  [meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

  mutil_c         =  par.nu*mutil_c_a(:)+(1-par.nu)*mutil_c_n;


distN=100;
countN=0;

while distN>mpar.crit && countN<100
    countN=countN+1;
    YminusNEW = (TFP*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));
    Rminus =  TFP *(1-par.alpha)   *mc.* ((Nminus)./Kminus).^(par.alpha)  - par.delta;
    Wminus =  TFP *par.alpha       *mc.* (Kminus./(Nminus)).^(1-par.alpha);

    Profitminus = (1-mc)*YminusNEW - YminusNEW.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2 + 1/2*par.phi*((K-Kminus).^2)./Kminus;

        % Incomes (grids)
    inc.rent    = meshes.k*Rminus;
    inc.capital = meshes.k*Qminus;
    inc.money   = (RBminus/PIminus).*meshes.m...
        + (meshes.m<0).*(par.borrwedge/PIminus).*meshes.m;
    
    EVk = reshape(reshape(Vk,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);
    RBaux = exp(RB)/PI + (meshes.m<0).*(par.borrwedge/PI);
    EVm = reshape(reshape(RBaux(:).*mutil_c,[mpar.nm*mpar.nk mpar.nh])*P',[mpar.nm,mpar.nk mpar.nh]);
    
    [c_a_star,m_a_star,k_a_star,c_n_star,m_n_star] = EGM_policyupdate(EVm,EVk,Qminus,PIminus,RBminus,Wminus,Profitminus,Tauminus,mc,Kminus,Nminus,inc,meshes,grid,par,mpar);
    
    naminus= ((1/par.gamma).* mutil(c_a_star(:)).*meshes.h(:).*par.tau.*par.alpha.*mc.*(Kminus./Nminus).^(1-par.alpha)).^(par.sigma_n);
    naminus          = min(naminus,1);
    naminus          = reshape(naminus,[mpar.nm mpar.nk mpar.nh]);
    naminus(:,:,end) = 0;
    nnminus= ((1/par.gamma).* mutil(c_n_star(:)).*meshes.h(:).*par.tau.*par.alpha.*mc.*(Kminus./Nminus).^(1-par.alpha)).^(par.sigma_n);
    nnminus          = min(nnminus,1);
    nnminus          = reshape(nnminus,[mpar.nm mpar.nk mpar.nh]);
    nnminus(:,:,end) = 0;
    
    Nguess  = JDminus(:)'*(par.nu.*naminus(:)+(1-par.nu)*nnminus(:));
    
    distN = abs(Nminus-Nguess);
    
    Nminus=0.25*Nguess+0.75*Nminus;
end


%% Update Value Functions

AProb=par.nu*ones([mpar.nm mpar.nk mpar.nh]);

%% Differences for distributions

% Transitions: find next smallest on-grid value for money and capital choices
[Dist_m_a,idm_a] = genweight(m_a_star,grid.m);
[Dist_m_n,idm_n] = genweight(m_n_star,grid.m);
[Dist_k,idk_a]   = genweight(k_a_star,grid.k);

% Iterate Distribution function forward
JD_new = JDIteration_mex(reshape(JDminus,[mpar.nm, mpar.nk, mpar.nh]),...
    int32(idm_n),int32(idm_a),int32(idk_a),Dist_m_n,Dist_m_a,Dist_k,...
    reshape(AProb,[mpar.nm, mpar.nk, mpar.nh]),P);

JD = reshape(JD_new(:),[mpar.nm,mpar.nk,mpar.nh]);

%% Third Set: Government Budget constraint
% Return on bonds (Taylor Rule)
taxrevenue =(1-Tauminus).*Wminus.*Nminus +(1-Tauminus).*Profitminus;
      
B = grid.B* exp(par.rho_B * log((Bminus)/(targets.B)) ...
    + par.rho_B * log(RBminus/par.RB)...
    - (par.rho_B+par.gamma_pi) * log(PIminus/par.PI) ...
    - par.gamma_T * log((taxrevenue)/(targets.T)));
% B= Bminus*RBminus/PIminus+par.G-taxrevenue;

marginal_b=squeeze(sum(sum(JD,2),3));
marginal_k=squeeze(sum(sum(JD,1),3));
BNEW  = sum(grid.m.*marginal_b');
KNEW  = sum(grid.k.*marginal_k);

QuantDiff=[ KNEW-K,BNEW-B];

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

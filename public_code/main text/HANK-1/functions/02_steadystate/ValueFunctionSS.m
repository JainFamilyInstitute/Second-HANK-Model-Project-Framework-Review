function [AProb,V,H,EU]=ValueFunctionSS(m_star,c_star,P_H,mpar,par,grid)
%% VALUEFNUCTION_SS calculates the value function given policies.
%

% Copyright (c) 2014-02-28
% Christian Bayer, Ralph L�tticke, Lien Pham-Dao, and Volker Tjaden
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


% Initialize matrices
weight11  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight21  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight22  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);

%% find next smallest on-grid value for money and capital choices
[Dist_m_a,idm_a] = genweight(m_a_star,grid.m);
[Dist_m_n,idm_n] = genweight(m_n_star,grid.m);
[Dist_k,idk_a]   = genweight(cap_a_star,grid.k);
idk_n = repmat(ones(mpar.nm,1)*(1:mpar.nk),[1 1 mpar.nh]); %This is the actual point on grid


%% Transition matrix adjustment case
idm_a=repmat(idm_a(:),[1 mpar.nh]);
idk_a=repmat(idk_a(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nk*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:),idk_a(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:),idk_a(:)+1,idh(:));
index21 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:)+1,idk_a(:),idh(:));
index22 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:)+1,idk_a(:)+1,idh(:));


for hh=1:mpar.nh

    %Corresponding weights
    weight21_aux =  Dist_m_a(:,:,hh).*(1-Dist_k(:,:,hh));
    weight11_aux = (1-Dist_m_a(:,:,hh)).*(1-Dist_k(:,:,hh));
    weight22_aux =  (Dist_m_a(:,:,hh)).*(Dist_k(:,:,hh));
    weight12_aux =  (1-Dist_m_a(:,:,hh)).*(Dist_k(:,:,hh));

    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P_H(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P_H(hh,:);
    weight21(:,:,hh)=weight21_aux(:)*P_H(hh,:);
    weight22(:,:,hh)=weight22_aux(:)*P_H(hh,:);
end
%Dimensions (mxk,h,h')
weight11=permute(weight11,[1 3 2]);
weight22=permute(weight22,[1 3 2]);
weight12=permute(weight12,[1 3 2]);
weight21=permute(weight21,[1 3 2]);
rowindex=repmat(1:mpar.nm*mpar.nk*mpar.nh,[1 4*mpar.nh]);

H_a=sparse(rowindex(:),[index11(:); index21(:); index12(:); index22(:)],...
    [weight11(:); weight21(:); weight12(:); weight22(:)],mpar.nm*mpar.nk*mpar.nh,mpar.nm*mpar.nk*mpar.nh); % mu'(h',k'), a without interest

%% Policy Transition Matrix for no-adjustment case
weight11  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight21  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
idm_n=repmat(idm_n(:),[1 mpar.nh]);
idk_n=repmat(idk_n(:),[1 mpar.nh]);
index11 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_n(:),idk_n(:),idh(:));
index21 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_n(:)+1,idk_n(:),idh(:));


for hh=1:mpar.nh

    %Corresponding weights
    weight21_aux =  Dist_m_n(:,:,hh);
    weight11_aux = (1-Dist_m_n(:,:,hh));

    weight21(:,:,hh)=weight21_aux(:)*P_H(hh,:);
    weight11(:,:,hh)=weight11_aux(:)*P_H(hh,:);

end
weight11=permute(weight11,[1 3 2]);
weight21=permute(weight21,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nk*mpar.nh,[1 2*mpar.nh]);

H_n=sparse(rowindex,[index11(:); index21(:)],...
    [weight11(:); weight21(:)],mpar.nm*mpar.nk*mpar.nh,mpar.nm*mpar.nk*mpar.nh); % mu'(h',k'), a without interest

%% Joint transition matrix and transitions

APD_a=sparse(1:length(AProb(:)),1:length(AProb(:)),AProb(:));
APD_n=sparse(1:length(AProb(:)),1:length(AProb(:)),1-AProb(:));
H=APD_a*H_a + APD_n*H_n;

util = @(c)(c.^(1-par.sigma))./(1-par.sigma);

u=((1-AProb).*util(c_n_star)+AProb.*util(c_a_star));

distV=1;

% Logit
AC =  par.probit2.*((1-AProb(:)./par.nu).*log(1-AProb(:)./par.nu) + AProb(:)./par.nu.*log(AProb(:)./par.nu))...
    + par.probit1.*AProb(:)./par.nu -  (par.probit1-par.probit2*log(1+exp(par.probit1/par.probit2)));
AC = AC.*par.nu;


V =(speye(length(AC)) -par.beta*H)\(u(:)-AC);
count=1;
while distV>mpar.crit && count<1000
    count=count+1;
    Vnew=u(:)+par.beta.*H*V - AC;
    distV=max(abs(Vnew-V));
    V=Vnew;
end
V_a =util(c_a_star(:)) +par.beta*H_a*V;
V_n =util(c_n_star(:)) +par.beta*H_n*V;
V_a =reshape(V_a,[mpar.nm,mpar.nk,mpar.nh]);
V_n =reshape(V_n,[mpar.nm,mpar.nk,mpar.nh]);
V   =reshape(V  ,[mpar.nm,mpar.nk,mpar.nh]);

DV = V_a-V_n;
% Logit
AProb = 1./(1+exp(-((DV)-par.probit1)./par.probit2));

AProb=par.nu*max(min(AProb,1-1e-6),1e-6);

AProb=reshape(AProb,[mpar.nm mpar.nk mpar.nh]);

if nargout==4
    EU =(speye(length(AC)) -par.beta*H)\(u(:));
    count=1;
    distV=100;

    while distV>mpar.crit && count<1000
        count=count+1;
        EUnew=u(:)+par.beta.*H*EU;
        distV=max(abs(EUnew-EU));
        EU=EUnew;
    end
end
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

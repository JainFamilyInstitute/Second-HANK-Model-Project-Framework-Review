function [Poly,InvCheb,Gamma]=createSparseBasis(grid,mpar,maxdim,Xss)
% This Function creates a sparse Basis, POLY, for polynomial approximation.
% It assumes that the grid points coincide with the Chebyshev nodes for a
% relevant transformation of the state space/ basis functions.
% It also provides the generalized inverse, INVCHEB, to estimate
% coefficients of the polynomial from function valuies using least squares.
%
% In addition it provides a matrix GAMMA such that GAMMA*x is a
% perturbation of the marginal distribution functions stored in XSS
%
% Authors: Christian Bayer and Ralph C. Luetticke, Bonn and London.
%
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


Poly=[];
maxlevel=max([mpar.nm mpar.nh]);

%% Values of Chebyshev Polynomials (1st type) at Chebyshev Nodes
% basem = (exp(exp(exp(exp(linspace(0,log(log(log(log(1 +1)+1)+1)+1),mpar.nm))-1)-1)-1)-1);
% Tm = cos(pi* (0:maxlevel-1)' * basem)';

Tm = cos(pi* (0:maxlevel-1)' * (linspace(0.5/mpar.nm/2,1-0.5/mpar.nm*2,mpar.nm)))';
% Tk = cos(pi* (0:maxlevel-1)' * (linspace(0.5/mpar.nk*2,1-0.5/mpar.nk/2,mpar.nk)))';
Th = cos(pi* (0:maxlevel-1)' * (linspace(0.5/(mpar.nh-1),1-0.5/(mpar.nh-1),(mpar.nh-1))))';


for j1=1:(length(grid.h)-1)
%     for j2=1:length(grid.k)
        for j3=1:length(grid.m)
            if j1+j3<maxdim
                [TT1, TT3] = ndgrid(Tm(:,j3), [Th(:,j1); 0]);
                Poly=[Poly TT1(:).*TT3(:)];
            end
        end
%     end
end

for j2=1:length(grid.m)
%     for j3=1:length(grid.k)
        if j2<maxdim-1
            [TT1, TT3] = ndgrid(Tm(:,j2), [zeros(length(grid.h)-1,1);1]);
            Poly=[Poly TT1(:).*TT3(:)];
        end
%     end
end


InvCheb=(Poly'*Poly)\Poly';

%% Mapping for Histogram

Gamma=zeros(mpar.nm*mpar.nh,mpar.nm*mpar.nh-1-mpar.nh);
for j=1:mpar.nm*mpar.nh-1-mpar.nh
    Gamma(1:mpar.nm*mpar.nh-mpar.nh,j)=-Xss(1:mpar.nm*mpar.nh-mpar.nh);
    Gamma(j,j)=1-Xss(j);
    Gamma(j,j)=Gamma(j,j) -sum(Gamma(1:mpar.nm*mpar.nh-mpar.nh,j));
end

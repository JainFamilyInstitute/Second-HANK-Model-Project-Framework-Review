function [c_guess, inc]=policyguess(meshes,WW,RBRB,par,mpar)
%policyguess returns autarky policy guesses (in the first period only):
% c_a_guess, c_n_guess, psi_guess
% as well as income matrices later on used in EGM:
% inc
%
% Consumption is compositite leisure and physical consumption (x_it) in the
% paper, therefore labor income is reduced by the fraction of leisure
% consumed.


% Copyright (c) 2014-02-28
% Christian Bayer, Ralph L�tticke, Lien Pham-Dao, and Volker Tjaden
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

inc.labor   = par.tau1.*WW.*meshes.h;
inc.money   = RBRB.*meshes.m;
inc.profits = 0;%lump sum profits

%% Initial policy guesses:  Autarky policies as guess

% Consumption guess
c_guess = inc.labor + max(inc.money,0) + inc.profits;


end

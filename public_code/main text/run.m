%
% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Transmission of Monetary Policy with Heterogeneity in Household Portfolios', 
% by Ralph Luetticke
% https://sites.google.com/site/ralphluetticke/
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

%% Initialize workspace and load directories
clc
clearvars 
close all
%%

restoredefaultpath
oldpath = path;

Computername='HANK'

global TayRule_case

TayRule_case = 3;
if TayRule_case == 1
    com_fig_name ='SS_BASELINE'
elseif TayRule_case == 2
    com_fig_name = 'EmpGap_H21';
else
   com_fig_name = 'AsyEmp_H21';
end

save_IRF_data = com_fig_name ;
%save_IRF_data = [com_fig_name '_Positive'];

starttime=clock;
set(0,'defaulttextinterpreter','latex')

%% Produce HANK-2 results

cd('HANK-2')

mainskript

mainskript_decomposition

cd ..

%% Produce HANK-1 results

cd('HANK-1')

mainskript

cd ..

%% Produce RANK results

% cd('RANK')
% 
% % dynare has to be in Matlab path
% dynare NKIM_RepAgent.mod
% 
% cd ..


%% Plotting
% all plots can be found in ../saves

aggrshock           = 'MP';

load('HANK-2/SS_BASELINE.mat')

load(['saves/IRF_' save_IRF_data '_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat'])

load(['saves/' save_IRF_data '_IRF_LIQUID_data.mat'])

% plot_all_irfs


plot_all_irfs_RANK_excluded

% close all


exit



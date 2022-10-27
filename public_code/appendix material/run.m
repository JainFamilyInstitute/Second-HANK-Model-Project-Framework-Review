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
clear 
close all

restoredefaultpath
oldpath = path;

Computername='HANK'

global TayRule_case 

TayRule_case = 2; % 2; %3;

if TayRule_case == 1
    com_fig_name = 'Original';
elseif TayRule_case == 2
    com_fig_name = 'EmpGap';
else
   com_fig_name = 'AsymEmp';
end

starttime=clock;
set(0,'defaulttextinterpreter','latex')


%% Produce results for HANK-2 baseline for inequality dimensions
cd('HANK-2')
save_IRF_data = [com_fig_name '_HANK-2'];
mainskript

cd ..

clc
clear 
close all

%% For Appendix Models 
% 
% %% Produce results for HANK-2 Real debt
% 
% cd('HANK-2_real')
% save_IRF_data = [com_fig_name '_Positive_HANK-2_real'];
% mainskript
% 
% cd ..
% 
% clc
% clear 
% close all
% %% Produce results for HANK-2 25% portfolio adjustment
% 
% cd('HANK-2_liquid')
% save_IRF_data = [com_fig_name '_Positive_HANK-2_liquid'];
% mainskript
% 
% cd ..
% 
% clc
% clear 
% close all
% %% Produce results for HANK-2 scarce liquidity
% 
% cd('HANK-2_scarce')
% save_IRF_data = [com_fig_name '_Positive_HANK-2_scarce'];
% mainskript
% 
% cd ..
% 
% clc
% clear 
% close all
% %% Produce results for HANK-2 GHH preferences
% 
% cd('HANK-2_ghh')
% save_IRF_data = [com_fig_name '_Positive_HANK-2_ghh'];
% mainskript
% 
% cd ..
% 
% clc
% clear 
% close all
% %% Produce results for HANK-2 lump-sum profits
% 
% cd('HANK-2_profits')
% save_IRF_data = [com_fig_name '_Positive_HANK-2_profits'];
% mainskript
% 
% cd ..
% 
% clc
% clear 
% close all
%% Produce results for HANK-1

cd('HANK-1')

mainskript

cd ..

clc
clear 
close all

%% 
% %% Plotting
% plotting
% 
% close all
function [grid]=makegrid(mpar,grid)

%% Construct basic grids
% This section defines the basic grids for the problem.
k_min = 0;
k_max = 40*grid.K;
grid.k = (exp(exp(exp(exp(linspace(0,log(log(log(log(k_max - k_min+1)+1)+1)+1),mpar.nk))-1)-1)-1)-1+k_min);  % set up quadruple exponential grid
% 
m_min = -0.5;
m_max = 20*grid.K;
grid.m = (exp(exp(exp(exp(linspace(0,log(log(log(log(m_max - m_min+1)+1)+1)+1),mpar.nm))-1)-1)-1)-1+m_min);
grid.m(abs(grid.m)==min(abs(grid.m)))=0;

% k_min = 0;
% k_max = 40*grid.K;
% grid.k = (exp(linspace(0,log(k_max - k_min+1),mpar.nk))-1+k_min);  % set up quadruple exponential grid
% 
% m_min = -0.75;
% m_max = 20*grid.K;
% grid.m = (exp(exp(linspace(0,log(log(m_max - m_min+1)+1),mpar.nm))-1)-1+m_min);  % set up quadruple exponential grid
% grid.m(abs(grid.m)==min(abs(grid.m)))=0;


function [grid]=makegrid(mpar)


%% Quadruble Log Grid
% m_min = 0;
% m_max = 100;
% grid.m = (exp(exp(exp(exp(linspace(0,log(log(log(log(m_max - m_min+1)+1)+1)+1),mpar.nm))-1)-1)-1)-1+m_min);
% grid.m(abs(grid.m)==min(abs(grid.m)))=0;

m_min = 0;
m_max = 2000;
grid.m = (exp(exp(exp(exp(linspace(0,log(log(log(log(m_max - m_min+1)+1)+1)+1),mpar.nm))-1)-1)-1)-1+m_min);
grid.m(abs(grid.m)==min(abs(grid.m)))=0;


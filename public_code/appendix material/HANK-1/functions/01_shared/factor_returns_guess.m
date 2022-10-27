function [ N,w,Profits_fc,WW,RBRB,R_fc,Y,mc ] = factor_returns_guess(meshes,K,par,mpar)
%factor_returns
  
%%  GHH Preferences    

%% with Y(N)
    


%       w = ((log(par.PI) - par.beta*par.EPI)/par.kappa + par.mu);
%       
%       N = (par.tau1.*w).^(1/par.gamma);   
%       
%       Y = par.H*N;
%       
%       Profits_fc = (1-w)*Y;
  
%% with Y(K,N)    
    mc =  par.mu - (par.beta * par.EPI - log(par.PI))/par.kappa;
    N    =  0.5;%(par.tau1.*par.alpha.*K^(1-par.alpha).*mc).^(1/(1-par.alpha+par.gamma))/2;   
    w = par.alpha .*mc.* (K./N).^(1-par.alpha);
        
    R_fc = (1-par.alpha) .*mc.* (N./K).^(par.alpha)- par.delta;

 
    Y = (N).^(par.alpha).*K.^(1-par.alpha);
    Profits_fc = (1-mc)*Y;
 %%
    NW=(N/par.H).*w;
    WW=NW*ones(mpar.nm,mpar.ns*mpar.nh); %Wages
    WW(:,end)=Profits_fc*par.profitshare;
    RBRB = (1+R_fc+(meshes.m<0)*par.borrwedge);

end


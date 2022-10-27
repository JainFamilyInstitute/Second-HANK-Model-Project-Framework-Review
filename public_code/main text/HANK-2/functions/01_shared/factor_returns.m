function [ R_fc,W_fc,Profits_fc,WWa,WWn,RR,RBRB,Y ] = factor_returns(meshes,grid,N,n_a_guess,n_n_guess,par,mpar)
%factor_returns
  
    mc =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;
        
%%  GHH Preferences    
    W_fc = par.alpha .*mc.* (grid.K./N).^(1-par.alpha);
    
    % Before tax return on capital
    R_fc = (1-par.alpha) .*mc.* (N./grid.K).^(par.alpha)- par.delta;
 
    Y = (N).^(par.alpha).*grid.K.^(1-par.alpha);
    Profits_fc = (1-mc)*Y - Y.*(1/(1-par.mu))./par.kappa./2 .*log(par.PI).^2;
 
    NWa=(n_a_guess./par.H).*W_fc;
    WWa=NWa; %Wages
    WWa(:,:,end)=Profits_fc*par.profitshare;
    NWn=(n_n_guess./par.H).*W_fc;
    WWn=NWn; %Wages
    WWn(:,:,end)=Profits_fc*par.profitshare;
    RR = R_fc; %Rental rates
    RBRB = par.RB/par.PI + (meshes.m<0).*(par.borrwedge/par.PI);


end


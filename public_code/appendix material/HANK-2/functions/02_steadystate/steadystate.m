function [n_a_guess,n_n_guess,c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr,R_fc,W_fc,Profits_fc,Output,grid]=steadystate(P_H,grid,meshes,mpar,par)
%     Prepare items for EGM
%     Layout of matrices:
%       Dimension 1: money m
%       Dimension 2: capital k
%       Dimension 3: stochastic human capital h

joint_distr = ones(mpar.nm,mpar.nk,mpar.nh)/(mpar.nh*mpar.nk*mpar.nm);

KL=0.75*grid.K;
KH=1.25*grid.K;
iter=0;
    
excessKhandle= @(K)excessK(K,joint_distr, grid,P_H,mpar,par,meshes,iter);
 
options=optimset('Display','off','TolFun',mpar.crit,'TolX',mpar.crit);
tic   
[Kcand, excess]=fsolve(excessKhandle,grid.K,options);
toc

%% Update

[excess,n_a_guess,n_n_guess,c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr, R_fc,W_fc,Profits_fc,Output,N]=excessK(Kcand,joint_distr, grid,P_H,mpar,par,meshes,iter);
disp(excess)
grid.N=N;
end
function [excess,n_a_guess,n_n_guess,c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr, R_fc,W_fc,Profits_fc,Output,N]=excessK(K,joint_distr, grid,P_H,mpar,par,meshes,iter)
grid.K=K;
distN=9999;
count=0;

[ N,R_fc,W_fc,Profits_fc,WWa,RR,RBRB,Output] = factor_returns_guess(meshes,grid,par,mpar);
WWn=WWa;
[c_a_guess,c_n_guess, psi_guess,inc]=policyguess(meshes,WWa,WWn,RR,RBRB,par,mpar);


while distN>mpar.crit && count<60 
    count=count+1;
    
    [~,~, ~,inc]=policyguess(meshes,WWa,WWn,RR,RBRB,par,mpar);
    
    % Solve Policies and Joint Distribution
    disp('Solving household problem by EGM')
    tic
    [c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,distPOL]=...
        policies_SS(c_a_guess,c_n_guess,psi_guess, W_fc, Profits_fc,grid, inc, RR,RBRB,P_H,mpar,par,meshes);
    mutil_c_n = 1./(c_n_guess.^par.xi); % marginal utility at consumption policy no adjustment
    mutil_c_a = 1./(c_a_guess.^par.xi); % marginal utility at consumption policy adjustment
    toc
    
    disp('Calc Joint Distr')
    tic
    [joint_distr,distJD]=JDiteration(m_n_star,m_a_star,cap_a_star,P_H,par,mpar,grid);
    joint_distr=reshape(joint_distr,[mpar.nm mpar.nk mpar.nh]);
    toc

    n_a_guess= ((1/par.gamma)*mutil_c_a.*meshes.h.*par.tau.*W_fc).^(par.sigma_n);
    n_a_guess=min(n_a_guess,1);
    n_a_guess(:,:,end)=0;
    n_n_guess= ((1/par.gamma)*mutil_c_n.*meshes.h.*par.tau.*W_fc).^(par.sigma_n);
    n_n_guess=min(n_n_guess,1);
    n_n_guess(:,:,end)=0;   
            
    Nguess=(par.nu.*joint_distr(:)'*n_a_guess(:)+(1-par.nu).*joint_distr(:)'*n_n_guess(:));
    distN=abs(Nguess-N);
    N=real(0.25*Nguess+0.75*N);
    
    [R_fc,W_fc,Profits_fc,WWa,WWn,RR,RBRB,Output] = factor_returns(meshes,grid,N,n_a_guess,n_n_guess,par,mpar);
    
    
end

AggregateCapitalDemand=sum(grid.k.*sum(sum(joint_distr,1),3));
excess=(grid.K-AggregateCapitalDemand);

end

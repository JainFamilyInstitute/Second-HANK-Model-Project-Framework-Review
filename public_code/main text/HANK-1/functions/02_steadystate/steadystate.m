function [c_guess,n_guess,m_star,joint_distr,R_fc,W_fc,Profits_fc,Output,N,grid,inc]=steadystate(P_H,grid,mpar,par)
%     Prepare items for EGM
%     Layout of matrices:
%       Dimension 1: money m and capital k
%       Dimension 2: stochastic human capital h


[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

% 1) Construct relevant returns and interpolation matrices from KS LoMs
joint_distr = ones(mpar.nm,mpar.nh)/(mpar.nh*mpar.nm);

[ N,W_fc,Profits_fc,WW,RBRB] = factor_returns_guess(meshes,grid.K,par,mpar);

[c_guess,inc]=policyguess(meshes,WW,RBRB,par,mpar);

% 3) Guess initial policies and determine labor income
KL=0.75*grid.K;
KH=1.25*grid.K;

[excessL]=excessK(KL,c_guess,joint_distr, grid,P_H,mpar,par,meshes);

[excessH]=excessK(KH,c_guess,joint_distr, grid,P_H,mpar,par,meshes);

if sign(excessL)==sign(excessH)
    disp('ERROR! Sign not diff')
end

%% Brent
fa=excessL; fb=excessH;
a=KL;
b=KH;

if fa*fb>0
    error('f(a) und f(b) sollten unterschiedliche Vorzeichen haben');
end

c=a; fc=fa;   %Zu Beginn ist c = a

c=a; fc=fa; d=b-a; e=d;

iter=0;
maxiter=1000;

while iter<maxiter
    iter=iter+1
    
    if fb*fc>0
        c=a; fc=fa; d=b-a; e=d;
    end
    
    if abs(fc)<abs(fb)
        a=b; b=c; c=a;
        fa=fb; fb=fc; fc=fa;
    end
    
    tol=2*eps*abs(b)+mpar.crit; m=(c-b)/2; %Toleranz
    
    if (abs(m)>tol) && (abs(fb)>0) %Verfahren muss noch durchgeführt werden
        
        if (abs(e)<tol) || (abs(fa)<=abs(fb))
            d=m; e=m;
        else
            s=fb/fa;
            if a==c
                p=2*m*s; q=1-s;
            else
                q=fa/fc; r=fb/fc;
                p=s*(2*m*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            end
            if p>0
                q=-q;
            else
                p=-p;
            end
            s=e; e=d;
            if ( 2*p<3*m*q-abs(tol*q) ) && (p<abs(s*q/2))
                d=p/q;
            else
                d=m; e=m;
            end
        end
        
        a=b; fa=fb;
        
        if abs(d)>tol
            b=b+d;
        else
            if m>0
                b=b+tol;
            else
                b=b-tol;
            end
        end
    else
        break;
    end
    
    [fb,c_guess,~,joint_distr]=excessK(b,c_guess,joint_distr, grid,P_H,mpar,par,meshes);
    
end

Kcand=b;
grid.K=b;

%% Update

[excess,c_guess,n_guess,m_star,joint_distr,R_fc,W_fc,Profits_fc,Output,N,inc]=excessK(Kcand,c_guess,joint_distr, grid,P_H,mpar,par,meshes);

grid.N=N;
end
function [excess,c_guess,n_guess,m_star,joint_distr,R_fc,W_fc,Profits_fc,Output,N,inc]=excessK(K,c_guess,joint_distr, grid,P_H,mpar,par,meshes)
distN=9999;

[ N,W_fc,Profits_fc,WW,RBRB,R_fc,Output,mc] = factor_returns_guess(meshes,K,par,mpar);
par.RB=R_fc+1;
count=0;
[c_guess,inc]=policyguess(meshes,WW,RBRB,par,mpar);

while distN>mpar.crit && count<150
    count=count+1;
    
    [~,inc]=policyguess(meshes,WW,RBRB,par,mpar);
    
    
    % 5) Solve Policies and Joint Distribution
    disp('Solving household problem by EGM')
    tic
    [c_guess,m_star,distPOL]=...
        policies_SS(c_guess, grid, inc,RBRB,P_H,mpar,par,W_fc,Profits_fc,meshes);
    mutil_c = 1./(c_guess.^par.sigma);
    
    toc
    
    disp(([distPOL]));
    
    disp('Calc Joint Distr')
    
    [joint_distr]=JDiteration(m_star,P_H,mpar,grid);
    joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);
     
    n_guess= ((1/par.gamma)*mutil_c.*meshes.h.*par.tau1.*W_fc).^(par.sigma_n);
    n_guess=min(n_guess,1);
    n_guess(:,end) = 0;
      Nguess=joint_distr(:)'*n_guess(:);
    
    
    distN=abs(Nguess-N);
    N=0.25*Nguess+0.75*N;
    
    [W_fc,Profits_fc,WW,RBRB,R_fc,Output,mc] = factor_returns(meshes,K,N,n_guess,par,mpar);
    par.RB=R_fc+1;

    
end

AggregateCapitalDemand=grid.m*sum(joint_distr,2)-grid.B;
excess=(K-AggregateCapitalDemand);

end




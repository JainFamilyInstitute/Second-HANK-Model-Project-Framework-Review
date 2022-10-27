function [c_new,m_star,distPOL] = policies_SS(c_guess, grid, inc,RBRB,P,mpar,par,w,Profits_fc,meshes)
% POLICIES solves for the household policies for consumption, money and
% capital holdings 

%% Initializations
util  = @(c)(c.^(1-par.sigma))./(1-par.sigma);
mutil = @(c)(1./(c.^par.sigma));
invutil = @(u)(((1-par.sigma).*u).^(1/(1-par.sigma)));
invmutil = @(mu)((1./mu).^(1/par.sigma));

%% Apply EGM to solve for optimal policies and marginal utilities
money_expense  = repmat(grid.m',[1 mpar.nh*mpar.ns]);
distC = 99999;
    options=optimset('Display','off','TolFun',mpar.crit,'TolX',mpar.crit);

count=0;
while max([distC])>mpar.crit && count<3000
    count=count+1;
    
    %% Step 1: Update policies for only money adjustment
    mutil_c = 1./(c_guess.^par.sigma); % marginal utility at consumption policy no adjustment
    
    mutil_c=RBRB.*mutil_c; %take return on money into account
    aux=reshape(permute(mutil_c,[2 1]),[mpar.ns*mpar.nh mpar.nm]);
    % form expectations
    EMU_aux = par.beta*permute(reshape(P*aux,[mpar.ns*mpar.nh mpar.nm]),[2 1]);
    
    c_aux = 1./(EMU_aux.^(1/par.sigma));
    
    n_guess= ((1/par.gamma)*mutil(c_aux).*meshes.h.*par.tau1.*w).^(par.sigma_n);
    n_guess=min(n_guess,1);
    n_guess(:,end) = 0;
    NW=(n_guess./par.H).*w;
    WW=NW.*ones(mpar.nm,mpar.ns*mpar.nh); %Wages
    WW(:,end)=Profits_fc*par.profitshare;
    inc.labor   = par.tau1.*WW.*meshes.h;
     
    % Take borrowing constraint into account
    [c_new,m_star]=EGM_Step1_b(grid,inc,money_expense,c_aux,mpar,par,Profits_fc,w,meshes,n_guess,mutil,options);
    
    m_star(m_star>grid.m(end)) = grid.m(end);
    
    %% Step 6: Check convergence of policies
    distC = max((abs(c_guess(:)-c_new(:))));
    
    % Update c policy guesses
    c_guess=c_new;
    
end
distPOL=[distC];
end

%% SUBFUNCTIONS

function [c_update,m_update]=EGM_Step1_b(grid,inc,money_expense,c_aux,mpar,par,Profits_fc,w,meshes,n_guess,mutil,options)
%%EGM_Step1_b computes the optimal consumption and corresponding optimal money
% holdings in case the capital stock cannot be adjusted by taking the budget
% constraint into account.
% c_update(m,k,s*h,M,K):    Update for consumption policy under no-adj.
% m_update(m,k,s*h,M,K):    Update for money policy under no-adj.

%% EGM: Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice

m_star = (c_aux + money_expense - inc.labor - inc.profits);
RR = (par.RB+(m_star<0)*par.borrwedge)./par.PI;
m_star = m_star./RR;

% Identify binding constraints
binding_constraints = money_expense < repmat(m_star(1,:),[mpar.nm 1 ]);


%% Next step: Interpolate w_guess and c_guess from new k-grids
% using c(s,h,k',K), k(s,h,k'K

m_star = reshape(m_star,[mpar.nm mpar.ns*mpar.nh]);
c_aux= reshape(c_aux,[mpar.nm mpar.ns*mpar.nh]);

%Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
[c_update, m_update]=egm1b_aux_mex(grid.m,m_star,c_aux);
c_update = reshape(c_update,[mpar.nm, mpar.ns*mpar.nh]);
m_update = reshape(m_update,[mpar.nm, mpar.ns*mpar.nh]);

bbaux=find(binding_constraints);
for j=1:sum(binding_constraints(:))
    bbind=bbaux(j);
    fcons= @(c_constr)cons_at_constraint(c_constr,meshes.h(bbind),inc.money(bbind)-grid.m(1),par,mutil,w);
    [c_constr, excess]=fzero(fcons,c_update(bbind),options);
    c_update(bbind)=c_constr;
end
m_update(binding_constraints) = min(grid.m);

m_update(m_update>grid.m(end)) = grid.m(end);

end
function [ex_cons,c_const]=cons_at_constraint(c_const,hh,money,par,mutil,w)

n_guess= ((1/par.gamma)*mutil(c_const).*hh.*par.tau1.*w).^(par.sigma_n);
n_guess=min(n_guess,1);
WW=(n_guess./par.H).*w;
labor   = par.tau1.*WW.*hh;
ex_cons=c_const-(labor  + money);

end


function [c_star,m_star]=EGM_policyupdate(EVm,Rminus,Qminus,BtoK,BtoKminus,RBminus,PIminus,Hminus,Wminus,Kminus,Nminus,mc,Profitminus,Tauminus,inc,mesh,grid,par,mpar)
mutil = @(c)(1./(c.^par.sigma));
    options=optimset('Display','off','TolFun',mpar.crit,'TolX',mpar.crit);

%% EGM Step 1:
EMU = par.beta*reshape(EVm,[mpar.nm,mpar.nh]);%reshape(reshape(Vm,[mpar.nm*mpar.nk,mpar.nh])*P',[mpar.nm,mpar.nk,mpar.nh]);
c_new = 1./(EMU.^(1/par.sigma));

    
nminus=((1/par.gamma).*mutil(c_new).*mesh.h.*Tauminus.*par.alpha*mc.* (Kminus./(Nminus)).^(1-par.alpha)).^(par.sigma_n);
nminus          = min(nminus,1);
nminus          = reshape(nminus,[mpar.nm mpar.nh]);
nminus(:,end)   = 0;
WW=(nminus./Hminus).*Wminus;
WW(:,end)=Profitminus*par.profitshare;
% Incomes (grids)
inc.labor   = Tauminus*WW.*(mesh.h);

% Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_star_n = (c_new + (Qminus+BtoK)*mesh.m - inc.labor);

m_star_n = m_star_n./(Rminus+Qminus+BtoKminus*RBminus/PIminus);%+(m_star_n<0)*par.borrwedge/PIminus);

% Identify binding constraints
binding_constraints = mesh.m < repmat(m_star_n(1,:),[mpar.nm 1]);



% Next step: Interpolate w_guess and c_guess from new k-grids
% using c(s,h,k'), k(s,h,k'K

m_star_n = reshape(m_star_n,[mpar.nm mpar.nh]);
c_n_aux  = reshape(c_new,[mpar.nm mpar.nh]);

% Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
% Check monotonicity of m_star_n
if max(sum(abs(diff(sign(diff(m_star_n))))))~=0
    warning('non monotone future liquid asset choice encountered')
end

[c_n_star, m_n_star]=egm1b_aux_mex(grid.m,m_star_n,c_n_aux);

c_star = reshape(c_n_star,[mpar.nm, mpar.nh]);
m_star = reshape(m_n_star,[mpar.nm, mpar.nh]);

% Consumption when drawing assets m' to zero: Eat all Resources
% Resource = inc.labor + inc.money;
% c_star(binding_constraints) = Resource(binding_constraints)-Qminus*grid.m(1);
bbaux=find(binding_constraints);
for j=1:sum(binding_constraints(:))
    bbind=bbaux(j);
    fcons= @(c_constr)cons_at_constraint(c_constr,mesh.h(bbind),inc.money(bbind)-grid.m(1),Tauminus,par,mutil,par.alpha*mc.* (Kminus./(Nminus)).^(1-par.alpha));
    [c_constr, excess]=fzero(fcons,c_star(bbind),options);
    c_star(bbind)=c_constr;
end

m_star(binding_constraints) = min(grid.m);
m_star(m_n_star>grid.m(end)) = grid.m(end);

end
function [ex_cons,c_const]=cons_at_constraint(c_const,hh,money,Tauminus,par,mutil,w)

n_guess= ((1/par.gamma)*mutil(c_const).*hh.*Tauminus.*w).^(par.sigma_n);
n_guess=min(n_guess,1);
WW=(n_guess./par.H).*w;
labor   = Tauminus.*WW.*hh;
ex_cons=c_const-(labor  + money);

end

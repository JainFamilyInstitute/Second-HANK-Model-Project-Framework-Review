tic
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));

Xss=[squeeze(sum(sum(joint_distr,2),3)); ... % marginal distribution liquid
    squeeze(sum(sum(joint_distr,1),3))'; ... % marginal distribution illiquid
    squeeze(sum(sum(joint_distr,2),1)); ... % marginal distribution productivity
    log(par.RB); 0];

Yss=[invmutil(mutil_c_a(:)); invmutil(mutil_c_n(:)); invmutil(Vk(:)); log(par.Q); log(par.PI); log(Output);...
    log(par.tau); log(par.W) ; log(par.R); log(par.PROFITS); log(par.N); log(targets.T);...
    log(grid.K); log(targets.B); log(targets.C); log(1-par.mu)];

[Poly,InvCheb,Gamma2] = createSparseBasis(grid,mpar,mpar.maxdim,Xss);
n1 = size(Poly); % used for controls
n2 = size(Gamma2); %used for distributions


% Produce matrices to reduce state-space
oc = length(Yss) - 3*n1(1);

InvGamma                          = sparse(3*n1(1)+n2(2)+2+oc,3*n1(2)+n2(2)+2+oc);

Gamma_state                       = sparse(Gamma2);
InvGamma(1:n2(1)+2,1:n2(1)+2)                 = eye(n2(1)+2);

Gamma_control                                 = sparse(3*n1(1)+oc,3*n1(2)+oc);
Gamma_control(1:n1(1),1:n1(2))                = Poly;
InvGamma(n2(2)+2+(1:n1(1)),n2(2)+2+(1:n1(2))) = InvCheb';

Gamma_control(n1(1)+(1:n1(1)),n1(2)+(1:n1(2)))                  = Poly;
InvGamma(n2(2)+n1(1)+2+(1:n1(1)),n2(2)+n1(2)+2+(1:n1(2)))       = InvCheb';

Gamma_control(3*n1(1)+(1:oc),3*n1(2)+(1:oc))                  = eye(oc);
InvGamma(n2(2)+3*n1(1)+2+(1:oc),n2(2)+3*n1(2)+2+(1:oc))       = eye(oc);

InvGamma                                      = InvGamma';

numstates   = n2(2)+2;
numcontrols = 3*n1(2)+oc;
State       = zeros(numstates,1);
State_m     = State;
Contr       = zeros(numcontrols,1);
Contr_m     = Contr;

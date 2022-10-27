tic
invutil = @(u)(((1-par.sigma).*u).^(1/(1-par.sigma)));
invmutil = @(mu)((1./mu).^(1/par.sigma));

Xss=[joint_distr(:);...
    log(par.RB); log(grid.B/grid.K); 0; 0];

naux=(n_guess(:,1:mpar.nh-1));
Yss=[(c_guess(:)); log(par.Q); log(par.PI); log(Output);...
    log(par.gamma1); log(par.W) ; log(par.R); log(par.PROFITS); log(par.N);...
    log(targets.T); log(grid.K); log(grid.B); log(targets.C); log(par.mu)];


[Poly,InvCheb,Gamma2] = createSparseBasis(grid,mpar,mpar.maxdim,Xss);
n1 = size(Yss); % used for controls
n2 = size(Gamma2); %used for distributions


% Produce matrices to reduce state-space
oc = length(Yss) - mpar.nm*mpar.nh;
os = length(Xss) - mpar.nm*mpar.nh;

InvGamma                          = sparse(1*n1(1)+n2(2)+4+oc,1*n1(2)+n2(2)+4+oc);

Gamma_state                       = sparse(Gamma2);
InvGamma(1:n2(1)+4,1:n2(1)+4)                 = eye(n2(1)+4);

Gamma_control                                 = sparse(1*n1(1)+oc,1*n1(2)+oc);
% 
% Gamma_control(1:n1(1),1:n1(2))                = Poly;
% InvGamma(n2(2)+3+(1:n1(1)),n2(2)+3+(1:n1(2))) = InvCheb';
% 
% % Gamma_control(n1(1)+(1:n1(1)),n1(2)+(1:n1(2)))                  = Poly;
% % InvGamma(n2(2)+n1(1)+2+(1:n1(1)),n2(2)+n1(2)+2+(1:n1(2)))       = InvCheb';
% % Gamma_control(2*n1(1)+(1:n1(1)),2*n1(2)+(1:n1(2)))              = Poly;
% % InvGamma(n2(2)+2*n1(1)+2+(1:n1(1)),n2(2)+2*n1(2)+2+(1:n1(2)))   = InvCheb';
% 
% Gamma_control(1*n1(1)+(1:oc),1*n1(2)+(1:oc))                  = eye(oc);
% InvGamma(n2(2)+1*n1(1)+3+(1:oc),n2(2)+1*n1(2)+3+(1:oc))       = eye(oc);
% 
% InvGamma                                      = InvGamma';

numstates   = n2(2)+4;
numcontrols = n1(1);
State       = zeros(numstates,1);
State_m     = State;
Contr       = zeros(numcontrols,1);
Contr_m     = Contr;

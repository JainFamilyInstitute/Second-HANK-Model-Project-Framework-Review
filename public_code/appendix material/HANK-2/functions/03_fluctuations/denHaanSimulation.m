%% DenHaan Errors

[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

mpar.ns=3;
% Use Rouwenhorst method to approximate stochastic volatility
[P_S, grid.s] = rouwen(par.rhoS, 0, par.sigmaS/sqrt(1-par.rhoS^2), mpar.ns);
grid.s=exp(grid.s');
PS_S=cumsum(P_S,2);

x0=zeros(numstates,1);
MX=[eye(length(x0));gx];
IRF_state_sparse=[];
Prices=[];
x=x0;

mpar.maxlag=16;
epsilon=0*randn(RandStream('mcg16807', 'Seed',20180621),mpar.maxlag,1);
epsilon(1) = par.sigmaS;

JDminus=joint_distr;
RBminus=log(par.RB);
Krealized_dH=grid.K;
Brealized_dH=grid.B;
Hrealized_dH=par.H;
ss(1)=0;
XssRED=[squeeze(sum(sum(joint_distr(1:end,:,:),2),3)); ... % marginal distribution liquid
    squeeze(sum(sum(joint_distr(:,1:end,:),1),3))'; ... % marginal distribution illiquid
    squeeze(sum(sum(joint_distr(:,:,1:end),2),1)); ... % marginal distribution productivity
    log(par.RB); 0];

PIGamma= pinv(Gamma2);
for t=1:mpar.maxlag
    ControlNOW=(gx*x)';
    xFOR=hx*x;
    ControlFOR=(gx*xFOR)';
    
    [JDminus,RBminus,Prices(t,:),PF,W,naminus,nnminus,R,c_a,m_a,k_a,c_n,m_n,P]  = denHaansys_Prices(JDminus,RBminus,ss(t),...
        ControlFOR',  ControlNOW',...
        Yss,Gamma_control,par,mpar,grid,P_H,aggrshock,oc,targets);
    
     [ginibonds(t),ginicapital(t),giniearnings(t),giniwealth(t),giniconsumption(t),giniincome(t)] =...
                diststatsGINIS(JDminus,grid,mpar,par.nu*c_a(:)+(1-par.nu)*c_n(:),R,RBminus,Prices(1),W,naminus,nnminus,PF,meshes,par);
    
    Krealized_dH(t+1)=grid.k*sum(sum(JDminus,1),3)';
    JDHaux=squeeze(sum(sum(JDminus,1),2));
    Hrealized_dH(t+1)=grid.h(1:end-1)*JDHaux(1:end-1);
    Brealized_dH(t+1)=grid.m*sum(sum(JDminus,2),3);
    % liquid assets
    aux_m = squeeze(sum(sum(JDminus,2),3));
    % illiquid assets
    aux_k =  squeeze(sum(sum(JDminus,1),3))' ;
    % human capital
    aux_h = squeeze(sum(sum(JDminus,1),2));
    Xaux=([aux_m(1:end);aux_k(1:end);aux_h(1:end)]-XssRED(1:(mpar.nm+mpar.nk+mpar.nh)));
    
    x(1:end-2)=PIGamma*Xaux;
    x(end-1)=(RBminus)-log(par.RB);
    ss(t+1) = ss(t)*par.rhoS + par.sigmaS*epsilon(t);
    x(end)=ss(t+1);
 
    
end

%% Figure 1
folder='figure1\';

figurename=[com_fig_name '_IRF_GINIC_theta2_' num2str(100*par.theta_pi) '_theta3_' num2str(100*par.rho_R)   '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(1:mpar.maxlag-1,100*100*(giniconsumption(2:end)-giniconsumption(1)),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
xticks([0 4 8 12 16])
xlim([0 16])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName', 'arial','FontSize',25);
printpdf(gcf,['..\saves\' folder figurename])

figurename=[com_fig_name '_IRF_GINIW_theta2_' num2str(100*par.theta_pi) '_theta3_' num2str(100*par.rho_R)   '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(1:mpar.maxlag-1,100*100*(giniwealth(2:end)-giniwealth(1)),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
xticks([0 4 8 12 16])
xlim([0 16])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName', 'arial','FontSize',25);
printpdf(gcf,['..\saves\' folder figurename])

figurename=[com_fig_name '_IRF_GINIEarnings_theta2_' num2str(100*par.theta_pi) '_theta3_' num2str(100*par.rho_R)   '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(1:mpar.maxlag-1,100*100*(giniearnings(2:end)-giniearnings(1)),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
xticks([0 4 8 12 16])
xlim([0 16])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName', 'arial','FontSize',25);
printpdf(gcf,['..\saves\' folder figurename])

figurename=[com_fig_name '_IRF_GINIIncome_theta2_' num2str(100*par.theta_pi) '_theta3_' num2str(100*par.rho_R)   '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(1:mpar.maxlag-1,100*100*(giniincome(2:end)-giniincome(1)),'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
xticks([0 4 8 12 16])
xlim([0 16])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName', 'arial','FontSize',25);
printpdf(gcf,['..\saves\' folder figurename])


IRF_GC = 100*100*(giniconsumption(2:end)-giniconsumption(1));
IRF_GW = 100*100*(giniwealth(2:end)-giniwealth(1));
IRF_GE = 100*100*(giniearnings(2:end)-giniearnings(1));
IRF_GI = 100*100*(giniincome(2:end)-giniincome(1));

cd ..
save(['saves/' save_IRF_data],"IRF_GI","IRF_GE","IRF_GW","IRF_GC");

cd('HANK-2')

% %% DenHaan Errors
% 
% [meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);
% 
% mpar.ns=3;
% % Use Rouwenhorst method to approximate stochastic volatility
% [P_S, grid.s] = rouwen(par.rhoS, 0, par.sigmaS/sqrt(1-par.rhoS^2), mpar.ns);
% grid.s=exp(grid.s');
% PS_S=cumsum(P_S,2);
% 
% x0=zeros(numstates,1);
% MX=[eye(length(x0));gx];
% IRF_state_sparse=[];
% Prices=[];
% x=x0;
% 
% mpar.maxlag=1000;
% epsilon=randn(RandStream('mcg16807', 'Seed',20180621),mpar.maxlag,1);
% pr_s=rand(RandStream('mcg16807', 'Seed',20180621),1,mpar.maxlag);
% 
% JDminus=joint_distr;
% RBminus=log(par.RB);
% Krealized_dH=grid.K;
% Brealized_dH=grid.B;
% Hrealized_dH=par.H;
% ss(1)=0;
% XssRED=[squeeze(sum(sum(joint_distr(1:end,:,:),2),3)); ... % marginal distribution liquid
%     squeeze(sum(sum(joint_distr(:,1:end,:),1),3))'; ... % marginal distribution illiquid
%     squeeze(sum(sum(joint_distr(:,:,1:end),2),1)); ... % marginal distribution productivity
%     log(par.RB); 0];
% 
% PIGamma= pinv(Gamma2);
% for t=1:mpar.maxlag
%     ControlNOW=(gx*x)';
%     xFOR=hx*x;
%     ControlFOR=(gx*xFOR)';
%     
%     [JDminus,RBminus,Prices(t,:),PF,W,naminus,nnminus,R,c_a,m_a,k_a,c_n,m_n,P]  = denHaansys_Prices(JDminus,RBminus,ss(t),...
%         ControlFOR',  ControlNOW',...
%         Yss,Gamma_control,par,mpar,grid,P_H,aggrshock,oc,targets);
%     
%      [ginibonds(t),ginicapital(t),giniearnings(t),giniwealth(t),giniconsumption(t),giniincome(t)] =...
%                 diststatsGINIS(JDminus,grid,mpar,par.nu*c_a(:)+(1-par.nu)*c_n(:),R,RBminus,Prices(1),W,naminus,nnminus,PF,meshes,par);
%     
%     Krealized_dH(t+1)=grid.k*sum(sum(JDminus,1),3)';
%     JDHaux=squeeze(sum(sum(JDminus,1),2));
%     Hrealized_dH(t+1)=grid.h(1:end-1)*JDHaux(1:end-1);
%     Brealized_dH(t+1)=grid.m*sum(sum(JDminus,2),3);
%     % liquid assets
%     aux_m = squeeze(sum(sum(JDminus,2),3));
%     % illiquid assets
%     aux_k =  squeeze(sum(sum(JDminus,1),3))' ;
%     % human capital
%     aux_h = squeeze(sum(sum(JDminus,1),2));
%     Xaux=([aux_m(1:end);aux_k(1:end);aux_h(1:end)]-XssRED(1:(mpar.nm+mpar.nk+mpar.nh)));
%     
%     x(1:end-2)=PIGamma*Xaux;
%     x(end-1)=(RBminus)-log(par.RB);
%         ss(t+1) = ss(t)*par.rhoS + par.sigmaS*epsilon(t);
%         x(end)=ss(t+1);
% 
%     
%     
% end
% 
% %% Linear Simulation
% 
% mpar.maxlag=t;
% 
% xlinear=x0;
% IRF_state_sparse=[];
% 
% for tt=1:mpar.maxlag
%     xlinear(end)= ss(tt);
%     IRF_state_sparse(:,tt)=(MX*xlinear)';
%     xlinear=hx*xlinear;
% end
% 
% %% Plotting
% IRF_distr  = Gamma_state*IRF_state_sparse(1:numstates-2,1:mpar.maxlag);
% clear Stateguess
% % preparation
% Stateguess(1,:)=       IRF_state_sparse(n2(2)+2,:);
% Stateguess(2,:)=grid.m*IRF_distr((1:mpar.nm),:)    + grid.B;
% Stateguess(3,:)=grid.k*IRF_distr(mpar.nm+(1:mpar.nk),:) + grid.K;
% Stateguess(4,:)=grid.h(1:end-1)*IRF_distr(mpar.nm+mpar.nk+(1:mpar.nh-1),:) + par.H;
% Stateguess(5,:)=       IRF_state_sparse(n2(2)+1,:) + log(par.RB);
% 
% PIlinear=1+IRF_state_sparse(end-oc+2,1:end);
% Qlinear=par.Q*(1+IRF_state_sparse(end-oc+1,1:end));
% PIrealized=Prices(:,2);
% Qrealized=Prices(:,1);
% logPI=log(PIrealized(:)./PIlinear(:));
% logQ=log(Qrealized(:)./Qlinear(:));
% 
% Klinear=Stateguess(3,:);
% Blinear=Stateguess(2,:);
% Hlinear=Stateguess(4,:);
% 
% logK=log(Krealized_dH(1:end-1)./Klinear);
% logB=log(Brealized_dH(1:end-1)./Blinear);
% logH=log(Hrealized_dH(1:end-1)./Hlinear);
% %% Table 1
% 
% upto=1000;
% denHaanTable(1,1)=mean(100*[abs(logQ(1:upto))]);
% denHaanTable(2,1)=max(100*[abs(logQ(1:upto))]);
% denHaanTable(1,2)=mean(100*[abs(logK(1:upto))]);
% denHaanTable(2,2)=max(100*[abs(logK(1:upto))]);
% denHaanTable(1,3)=mean(100*[abs(logPI(1:upto))]);
% denHaanTable(2,3)=max(100*[abs(logPI(1:upto))]);
% denHaanTable(1,4)=mean(100*[abs(logB(1:upto))]);
% denHaanTable(2,4)=max(100*[abs(logB(1:upto))]);
% 
%  save('..\saves\table1\table1_input', 'denHaanTable')
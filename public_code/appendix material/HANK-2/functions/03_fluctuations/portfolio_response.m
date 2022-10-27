
%%
[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

c_0=par.nu.*c_a_star_ss+(1-par.nu).*c_n_star_ss;
m_0=par.nu.*m_a_star_ss+(1-par.nu).*m_n_star_ss;
k_0=par.nu.*k_a_star_ss+(1-par.nu).*meshes.k;

[lambdaLI_0,lambdaIL_0,lambdaC_0,lambdaM_0,lambdaK_0] = diststats_liquidity_weights(JD_ss,m_0,k_0,c_0,grid,1,1,mpar);


c_2=par.nu.*c_a_star_a2+(1-par.nu).*c_n_star_a2;
m_2=par.nu.*m_a_star_a2+(1-par.nu).*m_n_star_a2;
k_2=par.nu.*k_a_star_a2+(1-par.nu).*meshes.k;
[lambdaLI_1,lambdaIL_1,lambdaC_1,lambdaM_1,lambdaK_1] = diststats_liquidity_weights(JD_a1,m_2,k_2,c_2,grid,Q(quart),1,mpar);

%% Plotting

IRF_PORT_CROSS=100*((lambdaLI_1(11:end)./lambdaIL_1(11:end))-(lambdaLI_0(11:end)./lambdaIL_0(11:end)));
figurename=['IRF_PortLiquidityPPGE_theta2_' num2str(100*par.theta_pi)  '_rhoB_' num2str(100*par.rho_B) '_quarter_' num2str(quart) '_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(11:100,IRF_PORT_CROSS,'b-','LineWidth',4.5)
ylabel('Percentage points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',25)
hold on
ylim([-5 3])
xlim([10 100])
plot(11:100,zeros(1,length(11:100)),'--black'); xlim([10 100])
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.15 0.175 0.8 0.8])
% printpdf(gcf, ['..\saves\' figurename])

IRF_LIQUID_CROSS=100*(lambdaLI_1(11:end)./lambdaLI_0(11:end)-1);
figurename=['IRF_MoneyGE_theta2_' num2str(100*par.theta_pi)  '_rhoB_' num2str(100*par.rho_B) '_quarter_' num2str(quart) '_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(11:100,IRF_LIQUID_CROSS,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(11:100,zeros(1,length(11:100)),'--black'); xlim([10 100])
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.15 0.175 0.8 0.8])
% printpdf(gcf, ['..\saves\' figurename])

IRF_ILLIQUID_CROSS=100*(lambdaIL_1(11:end)./lambdaIL_0(11:end)-1);
figurename=['IRF_CapitalGE_theta2_' num2str(100*par.theta_pi)  '_rhoB_' num2str(100*par.rho_B) '_quarter_' num2str(quart) '_' casename];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(11:100,IRF_ILLIQUID_CROSS,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(11:100,zeros(1,length(11:100)),'--black'); xlim([10 100])
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.175 0.175 0.775 0.8])
% printpdf(gcf, ['..\saves\' figurename])

close all

IRFname=['../saves/PortfolioIRF_' casename '_theta2_' num2str(100*par.theta_pi) '_theta3_' num2str(100*par.rho_R)  '_gamma2_' num2str(100*par.gamma_T) '_ABS_' num2str(par.ABS*100) '_' aggrshock];
save(IRFname,'IRF_ILLIQUID_CROSS','IRF_LIQUID_CROSS','IRF_PORT_CROSS');



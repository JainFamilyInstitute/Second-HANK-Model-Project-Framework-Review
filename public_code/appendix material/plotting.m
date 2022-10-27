for ff=1:5
    if ff==1
        %% Figure 9
        folder='figure9/';
        Name='realdebt'
        mpar.maxlag= 40;
        legend_baseline='Baseline'
        legend_1='Real Debt'
        %%
        load('saves/IRF_SS_BASELINE_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')
        
        IRF_Y_Baseline=IRF_Y;
        IRF_C_Baseline=IRF_C;
        IRF_I_Baseline=IRF_I;
        IRF_G_Baseline=IRF_G;
        IRF_K_Baseline=IRF_K;
        IRF_M_Baseline=IRF_M;
        IRF_N_Baseline=IRF_N;
        IRF_W_Baseline=IRF_W;
        IRF_PROFITS_Baseline=IRF_PROFITS;
        IRF_S_Baseline=IRF_S;
        IRF_PI_Baseline=IRF_PI;
        IRF_RB_Baseline=IRF_RB;
        IRF_Q_Baseline=IRF_Q;
        IRF_D_Baseline=IRF_D;
        IRF_LP_Baseline=IRF_LP;
        
        %%
        load('saves/IRF_SS_BASELINE_REALDEBT_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')
        IRF_Y_minus=IRF_Y;
        IRF_C_minus=IRF_C;
        IRF_I_minus=IRF_I;
        IRF_G_minus=IRF_G;
        IRF_K_minus=IRF_K;
        IRF_M_minus=IRF_M;
        IRF_N_minus=IRF_N;
        IRF_PROFITS_minus=IRF_PROFITS;
        IRF_W_minus=IRF_W;
        IRF_S_minus=IRF_S;
        IRF_PI_minus=IRF_PI;
        IRF_RB_minus=IRF_RB;
        IRF_Q_minus=IRF_Q;
        IRF_D_minus=IRF_D;
        IRF_LP_minus=IRF_LP;
        
    elseif ff==2
        %% Figure 8
                folder='figure8/';
        Name='RB0'
        mpar.maxlag= 40;
        legend_baseline='Baseline'
        legend_1='Scarce Liquidity'
        
        load('saves/IRF_SS_BASELINE_scarce_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')
        IRF_Y_minus=IRF_Y;
        IRF_C_minus=IRF_C;
        IRF_I_minus=IRF_I;
        IRF_G_minus=IRF_G;
        IRF_K_minus=IRF_K;
        IRF_M_minus=IRF_M;
        IRF_N_minus=IRF_N;
        IRF_PROFITS_minus=IRF_PROFITS;
        IRF_W_minus=IRF_W;
        IRF_S_minus=IRF_S;
        IRF_PI_minus=IRF_PI;
        IRF_RB_minus=IRF_RB;
        IRF_Q_minus=IRF_Q;
        IRF_D_minus=IRF_D;
        IRF_LP_minus=IRF_LP;
    elseif ff==3
        %% Figure 14
                folder='figure14/';
        Name='PROFIT'
        mpar.maxlag= 40;
        legend_baseline='Baseline'
        legend_1='Lump-sum profits'
        
        load('saves/IRF_SS_BASELINE_PROFITS_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')
        IRF_Y_minus=IRF_Y;
        IRF_C_minus=IRF_C;
        IRF_I_minus=IRF_I;
        IRF_G_minus=IRF_G;
        IRF_K_minus=IRF_K;
        IRF_M_minus=IRF_M;
        IRF_N_minus=IRF_N;
        IRF_PROFITS_minus=IRF_PROFITS;
        IRF_W_minus=IRF_W;
        IRF_S_minus=IRF_S;
        IRF_PI_minus=IRF_PI;
        IRF_RB_minus=IRF_RB;
        IRF_Q_minus=IRF_Q;
        IRF_D_minus=IRF_D;
        IRF_LP_minus=IRF_LP;
        
    elseif ff==4
        %% Figure 19
                folder='figure19/';
        Name='rhoB0_G'
        mpar.maxlag= 40;
        legend_baseline='Baseline'
        legend_1='Balanced budget'
        
        load('saves/IRF_SS_BASELINE_rhoB0_G_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')
        IRF_Y_minus=IRF_Y;
        IRF_C_minus=IRF_C;
        IRF_I_minus=IRF_I;
        IRF_G_minus=IRF_G;
        IRF_K_minus=IRF_K;
        IRF_M_minus=IRF_M;
        IRF_N_minus=IRF_N;
        IRF_PROFITS_minus=IRF_PROFITS;
        IRF_W_minus=IRF_W;
        IRF_S_minus=IRF_S;
        IRF_PI_minus=IRF_PI;
        IRF_RB_minus=IRF_RB;
        IRF_Q_minus=IRF_Q;
        IRF_D_minus=IRF_D;
        IRF_LP_minus=IRF_LP;
        
    elseif ff==5
        %% Figure 20
                folder='figure20/';
        Name='rhoB0_T'
        mpar.maxlag= 40;
        legend_baseline='Baseline'
        legend_1='Balanced budget'
        
        load('saves/IRF_SS_BASELINE_rhoB0_T_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')
        IRF_Y_minus=IRF_Y;
        IRF_C_minus=IRF_C;
        IRF_I_minus=IRF_I;
        IRF_G_minus=IRF_G;
        IRF_K_minus=IRF_K;
        IRF_M_minus=IRF_M;
        IRF_N_minus=IRF_N;
        IRF_PROFITS_minus=IRF_PROFITS;
        IRF_W_minus=IRF_W;
        IRF_S_minus=IRF_S;
        IRF_PI_minus=IRF_PI;
        IRF_RB_minus=IRF_RB;
        IRF_Q_minus=IRF_Q;
        IRF_D_minus=IRF_D;
        IRF_LP_minus=IRF_LP;
    end
    
    %% Plotting
    figurename=['IRF_Y_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_Y_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_Y_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    set(gca,'XTick',[0 4 8 12 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    hold on
    plot([0 16],[0 0],'--black')
    legend({legend_baseline,legend_1},'Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_C_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_C_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_C_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.235 0.175 0.725 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_I_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_I_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_I_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    %%
    figurename=['IRF_G_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_G_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_G_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    legend({legend_baseline,legend_1},'Location','SouthEast')
    set(gca,'Position',[0.25 0.175 0.7 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_TAU_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],-IRF_G_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],-IRF_G_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    legend({legend_baseline,legend_1},'Location','SouthEast')
    set(gca,'Position',[0.25 0.175 0.7 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_N_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_N_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_N_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_W_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_W_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_W_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_PF_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_PROFITS_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_PROFITS_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    
    figurename=['IRF_PI_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_PI_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_PI_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend({legend_baseline,legend_1},'Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_RB_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_RB_Baseline(2:17),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_RB_minus(2:17),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    
    %%
    figurename=['IRF_Q_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_Q_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_Q_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend({legend_baseline,legend_1},'Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_r_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_D_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_D_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_LP_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:15],IRF_LP_Baseline(1:15),'b-','LineWidth',4.5)
    hold on
    plot([1:1:15],IRF_LP_minus(1:15),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_K_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],[0 IRF_K_Baseline(1:15)],'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],[0 IRF_K_minus(1:15)],'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend({legend_baseline,legend_1},'Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.26 0.175 0.7 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_M_robustness_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_M_Baseline(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_M_minus(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
end

%% Plotting OneAsset vs TwoAsset
% Figure 11 and 12

load('saves/IRF_SS_BASELINE_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

for ff=1:2
    
    if ff==1
        
                folder='figure11/';

        
        load('saves/SS_BASELINE_PHI0_IRF_LIQUID_data.mat')
        
        Name=['baseline_PHI0'];
        
    else
        
                folder='figure12/';
        
        load('saves/SS_BASELINE_PHI18_IRF_LIQUID_data.mat')
        
        Name=['baseline_recal'];
        
    end
    %% Dynare Plotting
    signchange=1;
    figurename=['IRF_Y_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_Y(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_Y(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    set(gca,'XTick',[0 4 8 12 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    hold on
    plot([0 16],[0 0],'--black')
    legend('heterog. portfolios','represent. portfolio','Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    print('test','-dpng','-r300')
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_C_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_C(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_C(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.235 0.175 0.725 0.795])
    printpdf(gcf, ['saves/' folder figurename])
        
    
    figurename=['IRF_I_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_I(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_I(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    %%
    figurename=['IRF_G_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_G(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_G(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend('heterog. portfolios','represent. portfolio','Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.25 0.175 0.7 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    %
    figurename=['IRF_TAU_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],-IRF_G(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],-IRF_LIQUID_G(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend('heterog. portfolios','represent. portfolio','Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.25 0.175 0.7 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_N_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_N(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_N(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_W_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_W(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_W(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_PF_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_PROFITS(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_PROFITS(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_PI_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_PI(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_PI(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend('heterog. portfolios','represent. portfolio','Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_RB_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_RB(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_RB(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    
    %%
    figurename=['IRF_Q_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_Q(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_Q(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend('heterog. portfolios','represent. portfolio','Location','SouthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_r_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_D(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_D(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_LP_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:15],IRF_LP(1:15),'b-','LineWidth',4.5)
    hold on
    xlim([0 16])
    ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_K_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],[0 IRF_K(1:15)],'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],[0 IRF_LIQUID_K(1:15)],'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    plot([0 16],[0 0],'--black')
    legend('two-asset','one-asset','Location','NorthEast')
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.26 0.175 0.7 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_M_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    plot([1:1:16],IRF_M(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_LIQUID_M(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
    figurename=['IRF_MK_ONEvsTWOASSET_' Name];
    figure('Name',figurename,'NumberTitle','off','Position',[1 1 600 600])
    IRF_MK=100*log((M./K)/MK);%(100*(M./K-grid.B/grid.K));
    IRF_MKLIQUID=100*log((M_Liquid./K_Liquid)/(MK_Liquid));%(100*(M./K-grid.B/grid.K));
    plot([1:1:16],IRF_MK(1:16),'b-','LineWidth',4.5)
    hold on
    plot([1:1:16],IRF_MKLIQUID(1:16),'r--o','LineWidth',4.5)
    xlim([0 16])
    ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
    xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
    set(gca,'XTick',[0 4 8 12 16])
    hold on
    set(gca, 'FontName','arial','FontSize',25);
    set(gca,'Position',[0.215 0.175 0.75 0.795])
    printpdf(gcf, ['saves/' folder figurename])
    
end

%%
folder='figure7/';
load('saves/PortfolioIRF_SS_BASELINE_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

IRF_PORT_CROSS_baseline=IRF_PORT_CROSS;

load('saves/PortfolioIRF_SS_BASELINE_liquid_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

IRF_PORT_CROSS_liquid=IRF_PORT_CROSS;


figurename=['IRF_PortLiquidity_Comparison'];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 700 700])
plot(11:100,IRF_PORT_CROSS_baseline,'b-','LineWidth',4.5)
hold on
plot(11:100,IRF_PORT_CROSS_liquid,'r.','LineWidth',4.5)
plot(11:100,0.1+zeros(1,90),'k-.','LineWidth',4.5)
ylabel('Percentage points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Percentile of net wealth','Interpreter','none','FontName','arial','FontSize',25)
hold on
ylim([-5 3])
xlim([10 100])
legend('Baseline 12.5%','25% adj. prob.','100% adj. prob.')
plot(11:100,zeros(1,length(11:100)),'--black'); xlim([10 100])
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.15 0.175 0.8 0.8])
printpdf(gcf, ['saves/' folder figurename])

%%
folder='figure10/';
load('saves/PortfolioIRF_SS_BASELINE_PHI0_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

figurename=['IRF_PortLiquidityPPGE_SS_BASELINE_PHI0' ];
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
printpdf(gcf, ['saves/' folder figurename])

%%
folder='figure13/';
load('saves/PortfolioIRF_SS_BASELINE_PROFITS_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

figurename=['IRF_PortLiquidityPPGE_SS_BASELINE_PROFITS' ];
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
printpdf(gcf, ['saves/' folder figurename])

%%
folder='figure18/';
load('saves/PortfolioIRF_SS_BASELINE_rhoB0_T_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

figurename=['IRF_PortLiquidityPPGE_SS_BASELINE_rhoB0_T' ];
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
printpdf(gcf, ['saves/' folder figurename])

load('saves/PortfolioIRF_SS_BASELINE_rhoB0_G_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat')

figurename=['IRF_PortLiquidityPPGE_SS_BASELINE_rhoB0_G' ];
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
printpdf(gcf, ['saves/' folder figurename])

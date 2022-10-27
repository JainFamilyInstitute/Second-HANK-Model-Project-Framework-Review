clc;clearvars;close all


% First Figures are with 2-asset HANK (then with one-asset HANK)

%% load the data
cd("main text/saves/")

AggH2_Original = load("IRF_SS_BASELINE_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat");

AggH2_AsyEmp = load("IRF_AsyEmp_H21_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat");

AggH2_EmpGap = load("IRF_EmpGap_H21_theta2_150_theta3_80_gamma2_0_ABS_0_MP.mat");

% AggH1_Original = load("SS_BASELINE_IRF_LIQUID_data.mat");
% 
% AggH1_AsyEmp = load("AsyEmp_H21_IRF_LIQUID_data.mat");
% 
% AggH1_EmpGap = load("EmpGap_H21_IRF_LIQUID_data.mat");

cd ..
cd ..

cd("appendix material/saves/")


IneqH2_Original = load("OriginalHANK-2.mat");

IneqH2_AsymGap = load("AsymEmpHANK-2.mat");

IneqH2_EmpGap = load("EmpGapHANK-2.mat");

cd .. 
cd ..
%% Plot the figures - Aggregate for HANK 2
mkdir("Reqd_Figures/Aggregate")

cd("Reqd_Figures/Aggregate")


count = 0;

tpoints = 1:24;


count=count+1;
figure(count);
hold on;
plot(tpoints,AggH2_Original.IRF_C(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,AggH2_AsyEmp.IRF_C(tpoints), 'b--','LineWidth',1.5),grid;
% plot(tpoints,AggH2_EmpGap.IRF_C(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Consumption for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 4 8 12 16 20 24])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,AggH2_Original.IRF_N(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,AggH2_AsyEmp.IRF_N(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,AggH2_EmpGap.IRF_N(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Employment for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 4 8 12 16 20 24])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,AggH2_Original.IRF_N(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,AggH2_AsyEmp.IRF_N(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,AggH2_EmpGap.IRF_N(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Employment for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 4 8 12 16 20 24])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,AggH2_Original.IRF_N(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,AggH2_AsyEmp.IRF_N(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,AggH2_EmpGap.IRF_N(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Employment for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 4 8 12 16 20 24])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,AggH2_Original.IRF_Y(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,AggH2_AsyEmp.IRF_Y(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,AggH2_EmpGap.IRF_Y(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Output for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 4 8 12 16 20 24])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,AggH2_Original.IRF_I(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,AggH2_AsyEmp.IRF_I(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,AggH2_EmpGap.IRF_I(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Investment for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 4 8 12 16 20 24])
saveas(gcf,strcat(t.String,'.png'));

%% Plot the figures - Ginis for HANK 2
cd ..
cd ..


mkdir("Reqd_Figures/Inequality-Positive")

cd("Reqd_Figures/Inequality-Positive")

tpoints = 1:15;
count=count+1;
figure(count);
hold on;
plot(tpoints,IneqH2_Original.IRF_GW(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,IneqH2_AsymGap.IRF_GW(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,IneqH2_EmpGap.IRF_GW(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Wealth Gini for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 3 6 9 12 15])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,IneqH2_Original.IRF_GC(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,IneqH2_AsymGap.IRF_GC(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,IneqH2_EmpGap.IRF_GC(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Consumption Gini for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 3 6 9 12 15])
saveas(gcf,strcat(t.String,'.png'));


count=count+1;
figure(count);
hold on;
plot(tpoints,IneqH2_Original.IRF_GI(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,IneqH2_AsymGap.IRF_GI(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,IneqH2_EmpGap.IRF_GI(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Income Gini for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 3 6 9 12 15])
saveas(gcf,strcat(t.String,'.png'));

count=count+1;
figure(count);
hold on;
plot(tpoints,IneqH2_Original.IRF_GE(tpoints), 'm--','LineWidth',1.5),grid;
plot(tpoints,IneqH2_AsymGap.IRF_GE(tpoints), 'b--','LineWidth',1.5),grid;
%plot(tpoints,IneqH2_EmpGap.IRF_GE(tpoints), 'g:','LineWidth',1.5),grid;
plot(tpoints,zeros(size(tpoints)),'k-','LineWidth',1);
grid on;
legend({'Original' 'AsymEmpGap' },'Location','Best','Interpreter','latex');
xlim([1 tpoints(end)]);
% ylim([-1.5 0.5]);
hold off;
t=title('IRF For Earnings Gini for HANK with 2 assets');
t.FontSize=14;
ylabel('Deviation (\%)', 'interpreter','latex','FontSize',14);
xlabel('Quarters', 'interpreter','latex','FontSize',14);
set(gca,'FontSize',14);
set(gca,'XTick',[0 3 6 9 12 15])
saveas(gcf,strcat(t.String,'.png'));



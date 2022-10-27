%%
% Money
money_pdf=sum(sum(joint_distr,2),3);
money_cdf=cumsum(money_pdf);
median_money=sum(money_cdf<.5);
pdf_condmoney=squeeze(sum(joint_distr,1));

% Capital
capital_pdf=sum(sum(joint_distr,1),3)';
capital_cdf=cumsum(capital_pdf)';
median_capital=sum(capital_cdf<.5);
pdf_condcapital=squeeze(sum(joint_distr,2));

% Human Capital
humancapital_pdf=squeeze(sum(sum(joint_distr,1),2));
humancapital_cdf=squeeze(cumsum(humancapital_pdf));
median_humancapital=sum(humancapital_cdf<.5);
pdf_condhumancapital=squeeze(sum(joint_distr,3));

% Total wealth
for k = 1:mpar.nk
    for m = 1:mpar.nm
        mplusk(m+(k-1)*mpar.nm) = grid.m(m)+grid.k(k);
    end
end
[mplusk, IX]=sort(mplusk);
moneycapital_pdf=sum(joint_distr,3);
moneycapital_pdf=moneycapital_pdf(IX);
moneycapital_cdf=cumsum(moneycapital_pdf);

[meshes.m,meshes.k,meshes.h]=ndgrid(grid.m,grid.k,grid.h);

T_s=1;
NW=N(T_s).*W(T_s)/(par.N.*par.W)-1;
WWWW=repmat(NW,[mpar.nm mpar.nk mpar.nh]);
WWWW(:,:,end)=PROFITS(T_s)/par.PROFITS-1;
earnings=par.N*par.W*meshes.h/par.H;
earnings(:,:,end)=par.PROFITS*par.profitshare;
assetincome=meshes.m*(par.RB/par.PI-1)+meshes.k*par.R;
consumption=(par.nu*c_a_guess+(1-par.nu)*c_n_guess);

%%
for ii=1:4
    
    if ii==1
        dim=1;
        asset=meshes.m;
        returns=(par.RB./PI(T_s)-par.RB./par.PI);
        
    elseif ii==2       
        dim=2;
        asset=meshes.k;
        returns=R(T_s)-par.R;
        
    elseif ii==3
        dim=3;
        asset=earnings.*WWWW;
        returns=1;
        
    elseif ii==4               
        dim=4;
        asset=meshes.k;
        returns=Q(T_s)-par.Q;
        
    end
    
    
    
    for t=1:1
        for j=1:mpar.nm
            for i=1:mpar.nk
                RD_condh(j,i,t)=squeeze(joint_distr(j,i,:)./sum(joint_distr(j,i,:)))' * squeeze(asset(j,i,:)*returns(t)./consumption(j,i,:));
            end
        end
    end
    RD_condh=100*RD_condh;
    RD_condh=RD_condh(IX);
    RD_condh=reshape(RD_condh,[mpar.nm*mpar.nk t]);
    
    
    moneycapital_quintiles=[1 sum(moneycapital_cdf<.2);sum(moneycapital_cdf<.2)+1 sum(moneycapital_cdf<.4);...
        sum(moneycapital_cdf<.4)+1 sum(moneycapital_cdf<.6); sum(moneycapital_cdf<.6)+1 ...
        sum(moneycapital_cdf<.8); sum(moneycapital_cdf<.8)+1 mpar.nk*mpar.nm; ];
    percentilesstring=['1. ';'2. ';'3. ';'4. ';'5. '];
    
    
    for t=1:1
        for i=1:5
            if i==1
                RD_wealthdistr(dim,i,t) = (moneycapital_pdf(moneycapital_quintiles(i,1):moneycapital_quintiles(i,2))/(moneycapital_cdf(moneycapital_quintiles(i,2))))'* ...
                    RD_condh(moneycapital_quintiles(i,1):moneycapital_quintiles(i,2),t);
            else
                RD_wealthdistr(dim,i,t) = (moneycapital_pdf(moneycapital_quintiles(i,1):moneycapital_quintiles(i,2))/(moneycapital_cdf(moneycapital_quintiles(i,2))-moneycapital_cdf(moneycapital_quintiles(i,1)-1)))'* ...
                    RD_condh(moneycapital_quintiles(i,1):moneycapital_quintiles(i,2),t);
            end
        end
    end
    clear RD_condh
    
end

%% Latex Table
filename=['redistribution_table_' casename '.tex'];
FID = fopen(['../saves/' filename], 'w');
fprintf(FID, '\\begin{table} \\caption{Exposure to monetary shocks by wealth holdings} \\label{Tab:RD} \n');
fprintf(FID, '\\begin{center}\n');
fprintf(FID, '\\begin{threeparttable}\n');
fprintf(FID, '\\begin{tabular*}{\\textwidth}{l @{\\extracolsep{\\fill}} cccc}\n');
fprintf(FID, ' \\hline \\hline \n');
fprintf(FID, '\\\\ &\\multicolumn{3}{c}{Income gains/losses} & Capital gains/losses   \\\\ \\cline{2-4} \n');
fprintf(FID, '                       & Interest  &  Dividends  & Labor/Profit & on real assets    \\\\\n');
fprintf(FID, 'By wealth quintiles     & $\\Delta (R^B_{t-1}/\\pi_t$)   & $\\Delta r_t$ &  $\\Delta (W_tN_t+\\Pi_t$) & $\\Delta q_t$   \\\\ \\hline \n');
fprintf(FID, '\\\\ 1.   &    %8.2f & %8.2f & %8.2f & %8.2f \n', RD_wealthdistr(:,1));
fprintf(FID, '\\\\ 2.  &    %8.2f & %8.2f & %8.2f & %8.2f \n', RD_wealthdistr(:,2));
fprintf(FID, '\\\\ 3.  &    %8.2f & %8.2f & %8.2f & %8.2f \n', RD_wealthdistr(:,3));
fprintf(FID, '\\\\ 4.  &    %8.2f & %8.2f & %8.2f & %8.2f \n', RD_wealthdistr(:,4));
fprintf(FID, '\\\\ 5.  &    %8.2f & %8.2f & %8.2f & %8.2f \n', RD_wealthdistr(:,5));
fprintf(FID, '\\\\ \\hline \n');
fprintf(FID, '\\end{tabular*} \n');
fprintf(FID, '\\begin{small}\n');
fprintf(FID, '\\begin{tablenotes}\n');
fprintf(FID, '\\item \\textit{Notes:} Gains and losses in percent of within group consumption in period 0 to a one-standard deviation monetary policy shock, $\\epsilon^D=9$ basis points. Results are expressed in terms of steady-state consumption and averaged by using frequency weights from the steady-state wealth distribution. \n');
fprintf(FID, '\\end{tablenotes}\n');
fprintf(FID, '\\end{small}\n');
fprintf(FID, '\\end{threeparttable}\n');
fprintf(FID, '\\end{center}\n');
fprintf(FID, '\\end{table}\n');
fclose(FID);


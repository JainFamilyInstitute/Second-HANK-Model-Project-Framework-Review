function [lambdaLI,lambdaIL,lambdaC,lambdaM,lambdaK] = diststats_liquidity_weights(joint_distr,m_prime,k_prime,c_prime,grid,Q,PI,mpar)

joint_distr=reshape(joint_distr, [mpar.nm mpar.nk mpar.nh]);
% mpar.nh=mpar.nh-1;
% joint_distr=joint_distr(:,:,1:end-1)./sum(sum(sum(joint_distr(:,:,1:end-1))));
[MESH.m,MESH.k,MESH.h]=meshgrid(grid.m,grid.k,grid.h);

for j=1:mpar.nm
    for i=1:mpar.nk
%         inclabor(j,i)=squeeze(joint_distr(j,i,1:end)./sum(joint_distr(j,i,1:end)))' * squeeze(inc.labor(j,i,1:end));                
        mfull(j,i)=squeeze(joint_distr(j,i,1:end)./sum(joint_distr(j,i,1:end)))' * squeeze(m_prime(j,i,1:end));
        kfull(j,i)=squeeze(joint_distr(j,i,1:end)./sum(joint_distr(j,i,1:end)))' * squeeze(k_prime(j,i,1:end));
        cfull(j,i)=squeeze(joint_distr(j,i,1:end)./sum(joint_distr(j,i,1:end)))' * squeeze(c_prime(j,i,1:end));
                
        mgrid(j,i)=squeeze(joint_distr(j,i,1:end)./sum(joint_distr(j,i,1:end)))' * squeeze(MESH.m(j,i,1:end));
        kgrid(j,i)=squeeze(joint_distr(j,i,1:end)./sum(joint_distr(j,i,1:end)))' * squeeze(MESH.k(j,i,1:end));


    end
end

[mesh.k,mesh.m]=meshgrid(Q*grid.k,grid.m);

mplusk=mesh.m(:)+mesh.k(:);
[mplusk, IX]=sort(mplusk);

moverk=mesh.m(:)./mesh.k(:);
moverk(1:mpar.nk)=NaN;
moverk=moverk(IX);
mgrid=mgrid(IX);
kgrid=kgrid(IX);
mfull=mfull(IX);
kfull=kfull(IX);
cfull=cfull(IX);

moneycapital_pdf=max(sum(joint_distr,3),0);
moneycapital_pdf=moneycapital_pdf(IX);
moneycapital_cdf=cumsum(moneycapital_pdf);


% Select percentiles
P=100;
% Select smoothing parameter
h=0.1;

lambda=ones(P,1);
lambdaLI=ones(P,1);
lambdaIL=ones(P,1);
lambdaC=ones(P,1);
lambdaLIIL=ones(P,1);


data=[mplusk(:) mfull kfull moneycapital_cdf moneycapital_pdf*100 cfull mgrid kgrid];
data(data(:,5)==0,:)=[]; %drop combinations with zero prob    
        Y=(data(:,2))./(data(:,3)); % minus: obtain li/il ratio
        YLI=(data(:,2)); % Liquid assets
        YIL=(data(:,3)); % Illiquid assets 
        C=(data(:,6)); % Illiquid assets     
        AP=(data(:,7)); % Illiquid assets   
        M=data(:,7);
        K=data(:,8);


for j=1:P   
        
        dist=data(:,4)-j/P-0.5/P;
        weight=sqrt(normpdf(dist/h)).*sqrt(data(:,5));
        X=[weight dist.*weight];


        beta=regress(Y.*weight,X);
        lambda(j)=beta(1);
        beta=regress(YLI.*weight,X);
        lambdaLI(j)=(beta(1));
        beta=regress(YIL.*weight,X);
        lambdaIL(j)=(beta(1));
        beta=regress(C.*weight,X);
        lambdaC(j)=(beta(1));
        beta=regress(M.*weight,X);
        lambdaM(j)=(beta(1));
        beta=regress(K.*weight,X);
        lambdaK(j)=(beta(1));
end

end
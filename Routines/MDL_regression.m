function BEH_analysis_regression(dir,fin,model)

%load raw data
load(sprintf('BEHdata-%s.mat',fin),'BEHdata')

%load summary data
load(sprintf('BEHsummary-%s.mat',fin),'BEHsummary')

% load parameter fits
load(fullfile(dir.dir_fits,strcat('fit-', model, '-', fin, '.mat')),'results');

%load simulated variables
load(fullfile(dir.dir_fits ,strcat('sim-', model, '-',fin, '.mat')),'RPEs');

%% Regression: calculate regression and plot behavior

alldata = [];
for iS = 1:length(BEHdata)
    
    dec1 = BEHdata(iS).dec1';
    dec1 = dec1(:);
    dec1(1:100:end) = []; %remove first decision
    
    dec2 = BEHdata(iS).dec2';
    dec2 = dec2(:);
    dec2(1:100:end) = []; %remove first decision
    
    Qp1 = [diff(RPEs(iS).Q1_p1,[],2); diff(RPEs(iS).Q2_p1,[],2)];
    Qp1(100:100:end) = [];
    Qp2 = [diff(RPEs(iS).Q1_p2,[],2); diff(RPEs(iS).Q2_p2,[],2)];
    Qp2(100:100:end) = [];
    
    %add evidence per decision
    evidence = RPEs(iS).evidence_resp;
    evidence = evidence(:);
    evidence(100:100:end,:) = []; %remove final evidence
    
    %add credit assignment weight
    caw = RPEs(iS).cat;
    caw(1:100:end) = [];
    
    alldecs = [dec1; dec2];
    allcaw = [caw; caw];

    alldata = [alldata; alldecs Qp1 Qp2 evidence allcaw ones(594,1)*iS];
    
end
ds = dataset(alldata(:,1),zscore(alldata(:,2)),zscore(alldata(:,3)),zscore(alldata(:,4)),zscore(alldata(:,5)),alldata(:,6),'VarNames',{'choice','Qp1','Qp2','evidence','caw','VP'});

csvwrite(strcat('behaviorregression_',fin,'.csv'),double(ds))




%% prepare correlation/regression

wincrit = [50 50];
% wincrit = [60 40];
[fig1,mdl,ds] = calcswitch(fin,BEHdata,1,0,0,wincrit); %apply rescaling of rewards

nS = length(BEHdata);
[nB,nT] = size([BEHdata(1).dec1]);

rWiL = ismember(ds.Relev,1) & ismember(ds.Irrelev,-1); %relevantWin & irrelevantLoss
rLiW = ismember(ds.Irrelev,1) & ismember(ds.Relev,-1); %relevantLoss & irrelevantWin


for iS = 1:nS
    %subject mask
    vpident = ismember(ds.VP,iS);
    
    %get data for contrast
    dat_cond1 = double(ds(rWiL&vpident,1));
    dat_cond2 = double(ds(rLiW&vpident,1));
    
    %recode stay behavior
    dat1_cond1(dat_cond1==-1) = 0;
    dat1_cond2(dat_cond2==-1) = 0;
    
    %calculate average
    agg_cond1(iS) = sum(dat_cond1)/length(dat_cond1);
    agg_cond2(iS) = sum(dat_cond2)/length(dat_cond2);

end

% arcsine transformation of proportions
dall_cond1_trans = fw_asinetrans(agg_cond1);
dall_cond2_trans = fw_asinetrans(agg_cond2);


%% Correlation: correlate bandit performance, implicit CA and transfer performance

bandit = [BEHsummary.allCorr_block];
implicitCA = agg_cond1-agg_cond2;
transfer = [BEHsummary.allCorr_post];

x = zscore(results.x);
dsall = dataset(bandit',implicitCA',transfer',x(:,1),x(:,2),x(:,3),x(:,4),'VarNames',{'bandit','implicit','transfer','lr','beta','pers','ar'});

modelspec = 'implicit ~ lr + beta + pers + ar';
mdl = fitlm(dsall,modelspec)
reg_mean(:,1) = mdl.Coefficients.Estimate(2:end);
reg_sem(:,1) = mdl.Coefficients.SE(2:end);
reg_p(:,1) = mdl.Coefficients.pValue(2:end);

modelspec = 'bandit ~ lr + beta + pers + ar';
mdl = fitlm(dsall,modelspec)
reg_mean(:,2) = mdl.Coefficients.Estimate(2:end);
reg_sem(:,2) = mdl.Coefficients.SE(2:end);
reg_p(:,2) = mdl.Coefficients.pValue(2:end);

modelspec = 'transfer ~ lr + beta + pers + ar';
mdl = fitlm(dsall,modelspec)
reg_mean(:,3) = mdl.Coefficients.Estimate(2:end);
reg_sem(:,3) = mdl.Coefficients.SE(2:end);
reg_p(:,3) = mdl.Coefficients.pValue(2:end);



f2 = figure;
f2.Units = 'norm';
f2.Position = [0.5 0.1 0.5 0.3];

ax1 = axes;
% ax1.Position = [0.1 0.5 0.6 0.3];

hold on
b = bar([1 2 3 4],reg_mean,'grouped')

ngroups = 4;
nbars = 3;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    eb = errorbar(x, reg_mean(:,i), reg_sem(:,i), '.');
    eb.Color = 'k';
    eb.LineWidth = 2;
    
    
    
    for j = 1:ngroups
        txtY = max([reg_mean(j,i)+reg_sem(j,i)+0.02 0.02]);
        
        alphacorr = 1; %bonferroni = 12
       if reg_p(j,i)<.001/alphacorr
           text(x(j),txtY,'***','HorizontalAlignment','center','FontSize',16) %ax1.YLim(2)-0.02
       elseif reg_p(j,i)<.01/alphacorr
           text(x(j),txtY,'**','HorizontalAlignment','center','FontSize',16)
       elseif reg_p(j,i)<.05/alphacorr
           text(x(j),txtY,'*','HorizontalAlignment','center','FontSize',16)
       end
    end
end

leg = legend(b,{'Implicit credit assignment','Bandit performance','Transfer performance'})
leg.Location = 'eastoutside';
leg.Box = 'off';
title(leg,'Behavioral metric')

ax1.YLim = [-0.15 0.25];

xlabel('Model parameter')
ax1.XTick = [1 2 3 4];
ax1.XTickLabels = {'learning rate', 'inverse temperature', 'perseveration', 'assignment rate'}
ylabel('Standardized regression weight')

title(fin)
f2.Color = 'w';




%% psychometric plotting?
% nBin1 = 2;
% 
% ds2 = ds;
% mdl = fitglme(ds,'Qp1 ~ Qp2'); 
% ds2.Qp1 = mdl.Residuals.Raw;
% mdl = fitglme(ds,'Qp2 ~ Qp1'); 
% ds2.Qp2 = mdl.Residuals.Raw;
% 
% clear p1_rel p2_rel p1_rel_lowCAW p2_rel_lowCAW p1_rel_midCAW p2_rel_midCAW p1_rel_hihCAW p2_rel_hihCAW
% for iS = 1:length(BEHdata)
%     vpident = ismember(ds.VP,iS);
%     dsub = ds(vpident,:);
%     
% %     dsub2 = dsub;    
% %     mdl = fitglme(dsub,'Qp1 ~ Qp2'); 
% %     dsub2.Qp1 = mdl.Residuals.Raw;
% %     mdl = fitglme(dsub,'Qp2 ~ Qp1'); 
% %     dsub2.Qp2 = mdl.Residuals.Raw;    
% %     dsub = dsub2;
%     
%     %% for psychometric function
%     binlims = linspace(-2,2,nBin1+1);
%     Qp1bin = discretize(dsub.Qp1,binlims);
%     Qp2bin = discretize(dsub.Qp2,binlims);
%     
% %     quartiles = quantile(dsub.caw,[1/3 2/3]);
%     quartiles = quantile(dsub.caw,[1/2 1/2]);
%     idx1 = dsub.caw<=quartiles(1);
%     idx2 = dsub.caw>quartiles(1)&dsub.caw<quartiles(2);
%     idx3 = dsub.caw>=quartiles(2);
%     
%     for iBin = 1:nBin1
%         p1_rel(iS,iBin) = nanmean(dsub.choice(Qp1bin==iBin));
%         p2_rel(iS,iBin) = nanmean(dsub.choice(Qp2bin==iBin));
%         
%         p1_rel_lowCAW(iS,iBin) = nanmean(dsub.choice((Qp1bin==iBin)&idx1));
%         p2_rel_lowCAW(iS,iBin) = nanmean(dsub.choice((Qp2bin==iBin)&idx1));
%         
%         p1_rel_midCAW(iS,iBin) = nanmean(dsub.choice((Qp1bin==iBin)&idx2));
%         p2_rel_midCAW(iS,iBin) = nanmean(dsub.choice((Qp2bin==iBin)&idx2));
%         
%         p1_rel_hihCAW(iS,iBin) = nanmean(dsub.choice((Qp1bin==iBin)&idx3));
%         p2_rel_hihCAW(iS,iBin) = nanmean(dsub.choice((Qp2bin==iBin)&idx3));
%     end
%     
% end
% 
% %% psychometric plot
% x = binlims(1:end-1)+0.5.*diff(binlims);
% 
% f = figure;
% f.Units = 'norm';
% f.Position =[-0.9,0.3,0.9,0.5];
% 
% %% plot all dat
% ax = axes;
% ax.Position = [0.1 0.1 0.35 0.8];
% hold on
% 
% y = nanmean(p1_rel_lowCAW,1);
% [lerr] = fw_cousineau(p1_rel_lowCAW,'ci95');
% l1 = shadedErrorBar(x,y',lerr,'k',1);
% l1.mainLine.LineWidth = 2;
% 
% % y = nanmean(p1_rel_midCAW,1);
% % [lerr] = fw_cousineau(p1_rel_midCAW,'ci95');
% % l2 = shadedErrorBar(x,y',lerr,'b',1);
% % l2.mainLine.LineWidth = 2;
% 
% y = nanmean(p1_rel_hihCAW,1);
% [lerr] = fw_cousineau(p1_rel_hihCAW,'ci95');
% l3 = shadedErrorBar(x,y',lerr,'r',1);
% l3.mainLine.LineWidth = 2;
% 
% l4 = plot([binlims(2:end-1);binlims(2:end-1)],repmat([0; 1],1,size(binlims(2:end-1),2)),'--k','LineWidth',0.5);
% 
% % lgd = legend(ax,[l1.mainLine l2.mainLine],{'low' 'high'},'Location','southeast')
% % title(lgd,'arbitration')
% 
% ax.LineWidth = 2;
% ax.XTick = [0:20:100];
% ax.YTick = [0:.25:1];
% 
% xlabel('Q(A) minus Q(B)')
% ylabel('P(A)')
% 
% title('Correct policy')
% 
% %% arbitration weight - incorrect policy
% ax = axes;
% ax.Position = [0.5 0.1 0.35 0.8];
% hold on
% 
% y = nanmean(p2_rel_lowCAW,1);
% [lerr] = fw_cousineau(p2_rel_lowCAW,'ci95');
% l1 = shadedErrorBar(x,y',lerr,'k',1);
% l1.mainLine.LineWidth = 2;
% 
% % y = nanmean(p2_rel_midCAW,1);
% % [lerr] = fw_cousineau(p2_rel_midCAW,'ci95');
% % l2 = shadedErrorBar(x,y',lerr,'b',1);
% % l2.mainLine.LineWidth = 2;
% 
% y = nanmean(p2_rel_hihCAW,1);
% [lerr] = fw_cousineau(p2_rel_hihCAW,'ci95');
% l3 = shadedErrorBar(x,y',lerr,'r',1);
% l3.mainLine.LineWidth = 2;
% 
% l4 = plot([binlims(2:end-1);binlims(2:end-1)],repmat([0; 1],1,size(binlims(2:end-1),2)),'--k','LineWidth',0.5);
% 
% lgd = legend(ax,[l1.mainLine l3.mainLine],{'low'  'high'},'Location','southeast')
% % lgd = legend(ax,[l1.mainLine l2.mainLine l3.mainLine],{'low' 'medium' 'high'},'Location','southeast')
% title(lgd,'arbitration')
% lgd.Box = 'off';
% 
% ax.LineWidth = 2;
% ax.XTick = [0:20:100];
% ax.YTick = [0:.25:1];
% ax.YColor = 'none';
% 
% xlabel('Q(A) minus Q(B)')
% % ylabel('stay probability')
% 
% title('Incorrect policy')
% 


% 
% %% plot regression 
% 
% quartiles_evidence = quantile(ds.evidence,[1/3 2/3]);
% quartiles_caw = quantile(ds.caw,[1/3 2/3]);
% 
% modelspec = 'choice ~ Qp2 + Qp1 + (1 + Qp2 + Qp1 |VP)';
% 
% %negative evidence & high weight
% idx3 = ds.caw<=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% modelspec = 'choice ~ Qp1 + (1 + Qp1 |VP)';
% 
% %negative evidence & high weight
% idx3 = ds.caw<=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% modelspec = 'choice ~ Qp2 + (1 + Qp2 |VP)';
% 
% %negative evidence & high weight
% idx3 = ds.caw<=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% 
% modelspec = 'choice ~ Qp2 + Qp1 + (1 + Qp2 + Qp1 |VP)';
% 
% %negative evidence & high weight
% idx3 = ds.caw>=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% modelspec = 'choice ~ Qp1 + (1 + Qp1 |VP)';
% 
% %negative evidence & high weight
% idx3 = ds.caw>=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% modelspec = 'choice ~ Qp2 + (1 + Qp2 |VP)';
% 
% %negative evidence & high weight
% idx3 = ds.caw>=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% 
% 
% %% overall
% 
% ds2 = ds;
% 
% mdl = fitglme(ds,'Qp1 ~ Qp2'); 
% ds2.Qp1 = mdl.Residuals.Raw;
% 
% mdl = fitglme(ds,'Qp2 ~ Qp1'); 
% ds2.Qp2 = mdl.Residuals.Raw;
% 
% 
% nBin1 = 6;
% binlims = linspace(-2,2,nBin1+1);
% 
% Qp1bin = discretize(ds2.Qp1,binlims);
% Qp2bin = discretize(ds2.Qp2,binlims);
%     
% quartiles = quantile(ds2.caw,[1/3 2/3]);
% idx1 = ds2.caw<=quartiles(1);
% idx2 = ds2.caw>quartiles(1)&ds.caw<quartiles(2);
% idx3 = ds2.caw>=quartiles(2);
% 
% clear p1_rel p2_rel p1_rel_lowCAW p2_rel_lowCAW p1_rel_midCAW p2_rel_midCAW p1_rel_hihCAW p2_rel_hihCAW
% for iBin = 1:nBin1
%     p1_rel(iBin) = nanmean(ds2.choice(Qp1bin==iBin));
%     p2_rel(iBin) = nanmean(ds2.choice(Qp2bin==iBin));
%     
%     p1_rel_lowCAW(iBin) = nanmean(ds2.choice((Qp1bin==iBin)&idx1));
%     p2_rel_lowCAW(iBin) = nanmean(ds2.choice((Qp2bin==iBin)&idx1));
%     
%     p1_rel_midCAW(iBin) = nanmean(ds2.choice((Qp1bin==iBin)&idx2));
%     p2_rel_midCAW(iBin) = nanmean(ds2.choice((Qp2bin==iBin)&idx2));
%     
%     p1_rel_hihCAW(iBin) = nanmean(ds2.choice((Qp1bin==iBin)&idx3));
%     p2_rel_hihCAW(iBin) = nanmean(ds2.choice((Qp2bin==iBin)&idx3));
% end
% 
% %% psychometric plot
% x = binlims(1:end-1)+0.5.*diff(binlims);
% 
% f = figure;
% f.Units = 'norm';
% f.Position =[-0.9,0.3,0.9,0.5];
% 
% %% plot all dat
% ax = axes;
% ax.Position = [0.1 0.1 0.35 0.8];
% hold on
% 
% l1 = plot(x,p1_rel_lowCAW,'k','LineWidth',2);
% l3 = plot(x,p1_rel_hihCAW,'r','LineWidth',2);
% 
% 
% % y = nanmean(p1_rel_midCAW,1);
% % [lerr] = fw_cousineau(p1_rel_midCAW,'ci95');
% % l2 = shadedErrorBar(x,y',lerr,'b',1);
% % l2.mainLine.LineWidth = 2;
% 
% l4 = plot([binlims(2:end-1);binlims(2:end-1)],repmat([0; 1],1,size(binlims(2:end-1),2)),'--k','LineWidth',0.5);
% 
% % lgd = legend(ax,[l1.mainLine l2.mainLine],{'low' 'high'},'Location','southeast')
% % title(lgd,'arbitration')
% 
% ax.LineWidth = 2;
% ax.XTick = [0:20:100];
% ax.YTick = [0:.25:1];
% 
% xlabel('Q(A) minus Q(B)')
% ylabel('P(A)')
% 
% title('Correct policy')
% 
% %% arbitration weight - incorrect policy
% ax = axes;
% ax.Position = [0.5 0.1 0.35 0.8];
% hold on
% 
% l1 = plot(x,p2_rel_lowCAW,'k','LineWidth',2);
% l3 = plot(x,p2_rel_hihCAW,'r','LineWidth',2);
% 
% % y = nanmean(p2_rel_midCAW,1);
% % [lerr] = fw_cousineau(p2_rel_midCAW,'ci95');
% % l2 = shadedErrorBar(x,y',lerr,'b',1);
% % l2.mainLine.LineWidth = 2;
% 
% 
% l4 = plot([binlims(2:end-1);binlims(2:end-1)],repmat([0; 1],1,size(binlims(2:end-1),2)),'--k','LineWidth',0.5);
% 
% lgd = legend(ax,[l1 l3],{'low'  'high'},'Location','southeast')
% % lgd = legend(ax,[l1.mainLine l2.mainLine l3.mainLine],{'low' 'medium' 'high'},'Location','southeast')
% title(lgd,'arbitration')
% lgd.Box = 'off';
% 
% ax.LineWidth = 2;
% ax.XTick = [0:20:100];
% ax.YTick = [0:.25:1];
% ax.YColor = 'none';
% 
% xlabel('Q(A) minus Q(B)')
% % ylabel('stay probability')
% 
% title('Incorrect policy')
% 
% 
% 
% 
%     
% % modelspec = 'choice ~ Qp1 + Qp2 + (1 + Qp1 + Qp2 |VP)';
% % mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% % disp(mdl)
% % 
% % modelspec = 'choice ~ caw * (Qp1 + Qp2) + (1 + caw * (Qp1 + Qp2) |VP)';
% % mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% % disp(mdl)
% % 
% % modelspec = 'choice ~ evidence + (1 + evidence |VP)';
% % mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% % disp(mdl)
% % 
% % modelspec = 'choice ~ evidence * (Qp1 + Qp2) + (1 + evidence * (Qp1 + Qp2) |VP)';
% % mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% % disp(mdl)
% 
% fprintf('calculating full regression model...\n')
% modelspec = 'choice ~ evidence * (Qp1 + Qp2) + caw * (Qp1 + Qp2) + (1 + evidence * (Qp1 + Qp2) + caw * (Qp1 + Qp2)|VP)';
% modelspec = 'choice ~ evidence * (Qp1 + Qp2) + caw * (Qp1 + Qp2) + (1 |VP)';
% % modelspec = 'choice ~ evidence * caw * (Qp1 + Qp2) + (1 + evidence * caw * (Qp1 + Qp2)|VP)';
% mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl)
% 
% 
% %posthoc
% modelspec = 'choice ~ Qp1 + Qp2 + (1 + Qp1 + Qp2 |VP)';
% 
% %negative evidence
% quartiles = quantile(ds.evidence,[1/3 2/3]);
% idx1 = ds.evidence<=quartiles(1);
% dsub1 = ds(idx1,:);
% mdlhigh = fitglme(dsub1,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdlhigh)
% %positive evidence
% idx3 = ds.evidence>=quartiles(2);
% dsub2 = ds(idx3,:);
% mdllow = fitglme(dsub2,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdllow)
% 
% %low weights
% quartiles = quantile(ds.caw,[1/3 2/3]);
% idx1 = ds.caw<=quartiles(1);
% dsub1 = ds(idx1,:);
% mdlhigh = fitglme(dsub1,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdlhigh)
% %high weights
% idx3 = ds.caw>=quartiles(2);
% dsub2 = ds(idx3,:);
% mdllow = fitglme(dsub2,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdllow)
% 
% 
% %negative evidence & low weight
% quartiles_evidence = quantile(ds.evidence,[1/3 2/3]);
% quartiles_caw = quantile(ds.caw,[1/3 2/3]);
% idx1 = ds.evidence<=quartiles_evidence(1) & ds.caw<=quartiles_caw(1);
% dsub1 = ds(idx1,:);
% mdl_nElW = fitglme(dsub1,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nElW)
% 
% %positive evidence & low weight
% idx3 = ds.evidence>=quartiles_evidence(2) & ds.caw<=quartiles_caw(1);
% dsub2 = ds(idx3,:);
% mdl_pElW = fitglme(dsub2,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_pElW)
% 
% %negative evidence & high weight
% idx3 = ds.evidence<=quartiles_evidence(1) & ds.caw>=quartiles_caw(2);
% dsub3 = ds(idx3,:);
% mdl_nEhW = fitglme(dsub3,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_nEhW)
% 
% %positive evidence & high weight
% idx4 = ds.evidence>=quartiles_evidence(2) & ds.caw>=quartiles_caw(2);
% dsub4 = ds(idx4,:);
% mdl_pEhW = fitglme(dsub4,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
% disp(mdl_pEhW)
% 
% betaQs = [double(mdl_nElW.Coefficients(2:3,2))...
%     double(mdl_pElW.Coefficients(2:3,2))...
%     double(mdl_nEhW.Coefficients(2:3,2))...
%     double(mdl_pEhW.Coefficients(2:3,2))]'; 
% semQs = [double(mdl_nElW.Coefficients(2:3,3))...
%     double(mdl_pElW.Coefficients(2:3,3))...
%     double(mdl_nEhW.Coefficients(2:3,3))...
%     double(mdl_pEhW.Coefficients(2:3,3))]'; 
% 
% f = figure;
% f.Units = 'norm';
% f.Position = [0.1 0.3 0.3 0.4];
% f.Color = 'w';
% 
% ax1 = axes;
% ax1.Position = [0.1 0.15 0.4 0.7]
% ax1.LineWidth = 2;
% hold on
% b = bar([1 2],betaQs(1:2,:),'grouped')
% errorbar([0.85 1.85; 1.15 2.15]',betaQs(1:2,:),semQs(1:2,:),'LineStyle','none', 'Color', 'k','linewidth', 2)
% ax1.XTick = [1 2];
% ax1.XTickLabel = {'policy1' 'policy2'}
% ax1.YLim = [-0.50 7];
% title(ax1,'LOW')
% lgd = legend(ax1,b,{'V_{correct}' 'V_{incorrect}'})
% lgd.Location = 'north';
% lgd.FontSize = 10;
% lgd.Box = 'off';
% 
% ax2 = axes;
% ax2.Position = [0.55 0.15 0.4 0.7]
% ax2.LineWidth = 2;
% hold on
% b = bar([1 2],betaQs(3:4,:),'grouped')
% errorbar([0.85 1.85; 1.15 2.15]',betaQs(3:4,:),semQs(3:4,:),'LineStyle','none', 'Color', 'k','linewidth', 2)
% ax2.XTick = [1 2];
% ax2.XTickLabel = {'policy1' 'policy2'}
% ax2.YLim = [-0.50 7];
% ax2.YTick = [];
% ax2.YColor = 'none';
% title(ax2,'HIGH')
% 
% xlabel(ax1,'previous evidence')
% xlabel(ax2,'previous evidence')
% ylabel(ax1,'regression weight (a.u.)')
% 
% 
% txt1 = text(ax1,0,0,'arbitration weight')
% txt1.Position = [2.6,7.6,0];
% txt1.HorizontalAlignment = 'center';
% txt1.FontWeight = 'bold';
% 
% txt2 = text(2.75,-0.2,'arbitration weight = high')
% txt2.Position = [2.75,-0.2]
% txt2.HorizontalAlignment = 'center';










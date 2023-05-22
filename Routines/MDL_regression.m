function BEH_analysis_regression(dir,fin,model)

%load raw data
load(sprintf('BEHdata-%s.mat',fin),'BEHdata')

%load summary data
load(sprintf('BEHsummary-%s.mat',fin),'BEHsummary')

% load parameter fits
load(fullfile(dir.dir_fits,strcat('fit-', model, '-', fin, '.mat')),'results');

%% Regression: calculate regression and plot behavior
[~,~,ds] = calcswitch(fin,BEHdata,1,0,0); %apply rescaling of rewards

%% Contrast: ttest Stayprob(relWin & irLoss) vs Stayprob(relLoss & irWin)
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
       if reg_p(j,i)<.001
           text(x(j),ax1.YLim(2)-0.02,'***','HorizontalAlignment','center')
       elseif reg_p(j,i)<.01
           text(x(j),ax1.YLim(2)-0.02,'**','HorizontalAlignment','center')
       elseif reg_p(j,i)<.05
           text(x(j),ax1.YLim(2)-0.02,'*','HorizontalAlignment','center')
       end
    end
end

leg = legend(b,{'implicitCA','bandit','transfer'})
leg.Location = 'eastoutside';
leg.Box = 'off';
title(leg,'Behavioral metric')

ax1.YLim = [-0.08 0.2];

xlabel('Model parameter')
ax1.XTick = [1 2 3 4];
ax1.XTickLabels = {'learning rate', 'inverse temperature', 'perseveration', 'assignment rate'}
ylabel('Standardized regression weight')

title(fin)
f2.Color = 'w';







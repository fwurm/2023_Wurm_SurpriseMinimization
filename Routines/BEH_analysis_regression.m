function BEH_analysis_regression(exp_name)

%load raw data
load(sprintf('BEHdata-%s.mat',exp_name),'BEHdata')

%load summary data
load(sprintf('BEHsummary-%s.mat',exp_name),'BEHsummary')

%% Regression: calculate regression and plot behavior
wincrit = [40 60];
[fig1,mdl,ds] = calcswitch(exp_name,BEHdata,1,1,1,wincrit); %apply rescaling of rewards

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

% parametric test
fprintf('   Calculating paired ttest...\n')
[h_all,p_all,ci_all,stats_all] = ttest(dall_cond1_trans,dall_cond2_trans);
fprintf('      tstat = %.2f\n      pval = %.3f\n      df = %i\n',stats_all.tstat,p_all,stats_all.df)

% non-parametric test
fprintf('   Calculating paired Wilcoxon signed rank test...\n')
[p_all,h_all,stats_all] = signrank(dall_cond1_trans,dall_cond2_trans,'tail','both','method','approximate');
fprintf('      Z = %.2f\n      pval = %.3f\n',stats_all.zval,p_all)

%% Correlation: correlate bandit performance, implicit CA and transfer performance

bandit = [BEHsummary.allCorr_block];
implicitCA = agg_cond1-agg_cond2;
transfer = [BEHsummary.allCorr_post];

f2 = figure;
f2.Units = 'norm';
f2.Position = [-0.5 0.1 0.4 0.9];

ax1 = axes;
ax1.Position = [0.1 0.5 0.6 0.3];
scat = scatter(transfer,bandit); % correlation
scat.MarkerEdgeColor = 'k';
ls = lsline(ax1);
ls.LineWidth = 2;
ls.Color = 'k';
ax1.XLim = [0,1];
ax1.YLim = [0.4 1];
ax1.XTick = [0:0.2:1];
ax1.YTick = [0:0.2:1];
ax1.FontSize = 12;
ylabel('Bandit performance');

ax2 = axes;
ax2.Position = [0.1 0.1 0.6 0.3];
scat = scatter(transfer,implicitCA); % correlation
scat.MarkerEdgeColor = 'k';
ls = lsline(ax2);
ls.LineWidth = 2;
ls.Color = 'k';
ax2.XLim = [0,1];
ax2.YLim = [-0.5 1];
ax2.XTick = [0:0.2:1];
ax2.YTick = [-0.5:0.5:1];
ax2.FontSize = 12;
xlabel('Transfer performance');
ylabel('Implicit CA');

ax3 = axes;
ax3.Position = [0.1 0.82 0.6 0.1];
histo = histogram(transfer,linspace(0,1,13));
histo.EdgeColor = [0 0 0];
histo.FaceColor = [0.5 0.5 0.5];
ax3.XLim = [0 1];
ax3.XTick = [];
ax3.Box = 'off'
ax3.FontSize = 12;
ylabel('Frequency','FontSize',8);

ax4 = axes;
ax4.Position = [0.72 0.5 0.1 0.3];
histo = histogram(bandit,linspace(0.4,1,13));
histo.EdgeColor = [0 0 0];
histo.FaceColor = [0.5 0.5 0.5];
ax4.XLim = [0.4 1];
ax4.XTick = [];
ax4.Box = 'off';
ax4.FontSize = 12;
ax4.XDir = 'reverse';
% ylabel('Frequency','FontSize',8)
view([90 90]);

ax5 = axes;
ax5.Position = [0.72 0.1 0.1 0.3];
histo = histogram(implicitCA,linspace(-0.5,1,13));
histo.EdgeColor = [0 0 0];
histo.FaceColor = [0.5 0.5 0.5];
ax5.XLim = [-0.5 1];
ax5.XTick = [];
ax5.Box = 'off';
ax5.FontSize = 12;
ax5.XDir = 'reverse';
ylabel('Frequency','FontSize',8);
view([90 90]);

f2.Color = 'w'







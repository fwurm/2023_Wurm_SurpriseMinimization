function [f,mdl,ds] = calcswitch(model,BEHdata,rescale,doplot,doregress,wincrit)


f = [];
mdl = [];

if ~exist('rescale','var')
    rescale = 0;
end

if ~exist('doplot','var')
    doplot = 1;
end

if ~exist('doregress','var')
    doregress = 1;
end

if ~exist('wincrit','var')
    wincrit = [50 50];
end
wincritt = wincrit/50-1;

fprintf('######\n')
fprintf('## %s\n',model)
fprintf('######\n\n')

%number of subjects and blocks
nS = length(BEHdata);
[nB,nT] = size([BEHdata(1).dec1]);

%empty dataset to be appended per subject
data = [];
data2 = [];

% mark early and late trials
trialnums = repmat(1:nT-1,nB,1);
trialnum2 = reshape(trialnums',(nT-1)*nB,1);
trialstart = find(trialnum2<21);
trialend = find(trialnum2>=80);
trials = nan((nT-1)*nB,1);
trials(trialstart) = 1;
trials(trialend) = -1;

for iS = 1:nS
    
    % get stay behavior
    stays1 = double(~abs(diff(BEHdata(iS).dec1,[],2))); % Switch = 0; Stay = 1
    stays2 = double(~abs(diff(BEHdata(iS).dec2,[],2))); % Switch = 0; Stay = 1
    
    % get rewards
    R1 = squeeze(BEHdata(iS).rew1(:,1:end-1));
    R2 = squeeze(BEHdata(iS).rew2(:,1:end-1));
    
    if rescale
        R1 = (R1/50)-1;
        R2 = (R2/50)-1;
    end
    
    %code rewards as win or loss
    R1(R1<wincritt(1)) = -1;
    R1(R1>wincritt(2)) = 1;
    R1(R1>=wincritt(1)&R1<=wincritt(2)) = nan;
    
    R2(R2<wincritt(1)) = -1;
    R2(R2>wincritt(2)) = 1;
    R2(R2>=wincritt(1)&R2<=wincritt(2)) = nan;
    
    %reshape matrices to vectors
    stays1 = reshape(stays1',numel(stays1),1);
    stays2 = reshape(stays2',numel(stays1),1);
    R1 = reshape(R1',numel(stays1),1);
    R2 = reshape(R2',numel(stays1),1);
    
    %append data
    data = [data; stays1, R1, R2, ones(nB*(nT-1),1)*iS; stays2, R2, R1, ones(nB*(nT-1),1)*iS];
    data2 = [data2; stays1, R1, R2, ones(nB*(nT-1),1)*iS, trials; stays2, R2, R1, ones(nB*(nT-1),1)*iS, trials];
    
    losttrials(iS) = sum(sum(isnan(R1))) + sum(sum(isnan(R2)));
    
end

%create dataset
ds = dataset(data(:,1),data(:,2),data(:,3),data(:,4),'VarNames',{'Stays','Relev','Irrelev','VP'});
ds2 = dataset(data2(:,1),data2(:,2),data2(:,3),data2(:,4),data2(:,5),'VarNames',{'Stays','Relev','Irrelev','VP','startend'});

if doregress
    modelspec = 'Stays ~ Relev*Irrelev + (1 + Relev*Irrelev|VP)';
    mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
    disp(mdl)
    anova_all = anova(mdl)
    
    modelspec2 = 'Stays ~ Relev*Irrelev*startend + (1 + Relev*Irrelev*startend|VP)';
    mdl2 = fitglme(ds2,modelspec2,'Distribution','Binomial'); %Runs only on Matlab R2015b
    disp(mdl2)
    anova_all = anova(mdl2)
    
    startidx = ismember(ds2.startend,1);
    ds_start = ds2(startidx,:);
    mdl_start = fitglme(ds_start,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
    disp(mdl_start)
    anova_all = anova(mdl_start)
    
    endidx = ismember(ds2.startend,-1);
    ds_end = ds2(endidx,:);
    mdl_end = fitglme(ds_end,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
    disp(mdl_end)
    anova_all = anova(mdl_end)
end


% [B,BNAMES] = randomEffects(mdl);
% [h,p,ci,stats] = ttest(B(2:4:end),B(3:4:end))

if doplot
    vars = grpstats(ds,{'Relev' 'Irrelev'},{'mean','sem'});
    order = [4 3 2 1];
    
    f = figure;
    f.Units = 'norm';
    f.Position = [0.1 0.1 0.3 0.3];
    ax = axes;
    hold on
    % b = bar([1 2],vars.mean_Stays([4 3; 2 1]),'grouped'); %,'BaseValue',0.5
    b1 = bar([1 4],vars.mean_Stays([4 2])); %,'BaseValue',0.5
    b2 = bar([2 5],vars.mean_Stays([3 1])); %,'BaseValue',0.5
    
    
    b1.FaceColor = [0.3 0.3 0.3];
    b2.FaceColor = [0.8 0.8 0.8];
    b1.LineWidth = 2;
    b2.LineWidth = 2;
    
    b1.BarWidth = 0.3;
    b2.BarWidth = 0.3;
    
    a = errorbar([1 2 4 5],vars.mean_Stays([order]),vars.sem_Stays([order]),'.k','LineWidth',1);
    
    xlabel('Relevant outcome');
    ax.XLim = [0 6];
    ax.YLim = [0 1];
    
    ax.XTick = [1.5 4.5];
    ax.XTickLabels = {'win' 'loss'};
    
    leg = legend([b1,b2],{'win' 'loss'});
    title(leg,sprintf('Irrelevant\noutcome'));
    leg.Box = 'off';
    
    title(model);
    f.Color = 'w';
end

function [f,mdl,ds] = calcswitch(model,BEHdata,rescale,doplot,doregress)

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

fprintf('######\n')
fprintf('## %s\n',model)
fprintf('######\n\n')

%number of subjects and blocks
nS = length(BEHdata);
[nB,nT] = size([BEHdata(1).dec1]);

%empty dataset to be appended per subject
data = [];

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
    R1(R1<0) = -1;
    R1(R1>=0) = 1;
    R2(R2<0) = -1;
    R2(R2>=0) = 1;
    
    %reshape matrices to vectors
    stays1 = reshape(stays1',numel(stays1),1);
    stays2 = reshape(stays2',numel(stays1),1);
    R1 = reshape(R1',numel(stays1),1);
    R2 = reshape(R2',numel(stays1),1);
    
    %append data
    data = [data; stays1, R1, R2, ones(nB*(nT-1),1)*iS; stays2, R2, R1, ones(nB*(nT-1),1)*iS];
    
end

%create dataset
ds = dataset(data(:,1),data(:,2),data(:,3),data(:,4),'VarNames',{'Stays','Relev','Irrelev','VP'});

if doregress
    modelspec = 'Stays ~ Relev*Irrelev + (1 + Relev*Irrelev|VP)';
    mdl = fitglme(ds,modelspec,'Distribution','Binomial'); %Runs only on Matlab R2015b
    disp(mdl)
    anova_all = anova(mdl)
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

function batch_SurpriseMinimization
% in order to replicate the results from Wurm, Ernst, Steinhauser (2023), 
% download the preprocessed EEG data  and regression results via figshare 
% (https://figshare.com/articles/dataset/2023_Wurm_SurpriseMinimization/23055974).
% Unzip and change dir_data accordingly.

% This scripts requires following toolboxes to be installed (added to
% working directory):
%   - mfit (https://github.com/sjgershm/mfit.git)
%   - EEGLAB (https://github.com/sccn/eeglab.git)
%   - Mass Univariate ERP (https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox.git)

dir.dir = [pwd '\'];
dir.dir_model = [dir.dir 'Models\'];
dir.dir_fits = [dir.dir 'Fits\'];
dir.dir_simulation = [dir.dir 'Simulation\'];
dir.dir_eeg = 'C:\Users\wurmf\Documents\Data\eegCAT\EEG\';
% dir.dir_eeg = 'D:\Experiment_data\eegCAT_random\EEG\';
dir.setting = 'github'; %adapted for github

%adding model
addpath('Routines\');
addpath('Models\');
addpath('mfit\'); %add Gershman's mfit toolbox

%% Simulation

% models = {'hierarchicalRL-ar-counter'};
% models = {'hierarchicalRL-counter-corr90' 'hierarchicalRL-counter-incorr90'};
% models = {'jointRL'};

% SIM_simulateBehavior_rev(models,'empprior','value',1000) %model,priors,feedback,simulation

%% Behavioral data
% BEH_analysis_regression_rev('exp1-beh') %Fig 2ACE
% BEH_analysis_regression_rev('exp2-eeg') %Fig 2BDF

%% Model fitting
optionfit.type = 'behavioral fit';
optionfit.nstarts = [10];
optionfit.rewardscale = [0 100];
optionfit.nBlock = 3;
optionfit.nTrial_learn = 100;
optionfit.nTrial_post = 20;

% models = {'hierarchicalRL-ar' 'hierarchicalRL-ar-start' ...
%     'hierarchicalRL-ar-counter' 'hierarchicalRL-ar-counter-start' ...
%      'hierarchicalRL-ar-decay' 'hierarchicalRL-ar-decay-start' ...
%     'hierarchicalRL-corr90' 'hierarchicalRL-incorr90' ...
%     'hierarchicalRL-corr100' 'hierarchicalRL-incorr100' ...
%     'hierarchicalRL-rand' ...
%     'hierarchicalRL-counter-corr90' 'hierarchicalRL-counter-incorr90' ...
%     'hierarchicalRL-counter-corr100' 'hierarchicalRL-counter-incorr100' ...
%     'hierarchicalRL-counter-rand' ...
%     'jointRL'}; 

% for iM = 1:length(models) 
%     MDL_parafit(dir,'exp1-beh',models{iM},optionfit)
%     MDL_parafit(dir,'exp2-eeg',models{iM},optionfit)
% end

%% Model comparison

% model comparison in main text
models = {'hierarchicalRL-ar-counter',...
    'hierarchicalRL-counter-corr90',...
    'hierarchicalRL-counter-incorr90',...
    'hierarchicalRL-counter-rand',...
    'jointRL'};

% model comparison in supplement
% models = {'hierarchicalRL-ar', 'hierarchicalRL-ar-start',...
%     'hierarchicalRL-ar-decay', 'hierarchicalRL-ar-decay-start',...
%     'hierarchicalRL-ar-counter', 'hierarchicalRL-ar-counter-start'};

% model comparison for rebuttal
% models = {'hierarchicalRL-corr100' 'hierarchicalRL-incorr100',...
%     'hierarchicalRL-corr90' 'hierarchicalRL-incorr90',...
%     'hierarchicalRL-rand',...
%     'hierarchicalRL-counter-corr100' 'hierarchicalRL-counter-incorr100',...
%     'hierarchicalRL-counter-corr90' 'hierarchicalRL-counter-incorr90',...
%     'hierarchicalRL-counter-rand'};

% MDL_modelcomp(dir,'exp1-beh',models)
% MDL_modelcomp(dir,'exp2-eeg',models)



%% Simulate variables from fits
optionsim.type = 'simulation';
optionsim.nSamples = 10;
optionsim.fbtype = 'relpay';
optionsim.nSet = 2;
optionsim.nBlock = 3;
optionsim.nTrial_learn = 100; %number of simulated trials

model = 'hierarchicalRL-ar-counter';

% MDL_varsim(dir,'exp1-beh',model,optionsim)
% MDL_varsim(dir,'exp2-eeg',model,optionsim)


%% correlate behavior and parameters (new analysis based on review)

model = 'hierarchicalRL-ar-counter';

% MDL_regression(dir,'exp1-beh',model) %Fig 3A
% MDL_regression(dir,'exp2-eeg',model) %Fig 3B


%% single-trial regression FBLOCKED

EEGinfix = 'fbLocked\richef'; %EEG data for regression
outfix = 'fbLocked\regression'; %folder to save to [outfix model]

%settings for permutation test
timewin = [0 1000]; 
alpha = 0.05;
nPerm = 10000;


models = {'hierarchicalRL-ar-counter' 'hierarchicalRL-ar'  }; %'hierarchicalRL-ar-counter' 'hierarchicalRL-ar-decay' 'hierarchicalRL-ar-start'
regtypes = {'fb-revision' 'fb-revision-posthoc' 'fb-startend'}; 
for iM = models
    GNDinfix = strcat(outfix, '-', iM{:});    
    for i = 1:length(regtypes)
        EEG_regression(dir,'exp2-eeg',iM{:},regtypes{i},'in',EEGinfix,'out',outfix,'channels',[1:64],'srate',128,'parallel',4); %
        EEG_buildGND(dir,'exp2-eeg',iM{:},GNDinfix,regtypes{i},EEGinfix,[-0.5 1])
        EEG_permuteGND(dir,iM{:},outfix,regtypes{i},EEGinfix,timewin,alpha,nPerm)
    end
end

%plot
if 1
    model = 'hierarchicalRL-ar-counter';
    outfix = 'fbLocked\regression'; %folder to save to [outfix model]
    regtype = 'fb-revision-final'; %'fb-revision-final'
    GNDinfix = strcat(outfix, '-', model);
    
    %load results GND
    fin = fullfile(dir.dir_eeg, GNDinfix, regtype,['GND-' regtype '-perm.GND']);
    load(fin,'-mat');
    
    %plot to get an overview of effects
    for i = 1:length(GND_perm.t_tests)
        fprintf('##########\n')
        fprintf('# %s\n',GND_perm.bin_info(i).bindesc)
        fprintf('##########\n')
        fprintf('regression number: %d\n',i)
        fprintf('   positive clusters: %d\n',length(GND_perm.t_tests(i).clust_info.pos_clust_pval))
        fprintf('     significant: %d\n',length(find(GND_perm.t_tests(i).clust_info.pos_clust_pval<.05)))
        fprintf('     smallest pval: %.3f\n',min(GND_perm.t_tests(i).clust_info.pos_clust_pval))
        fprintf('   negative clusters: %d\n',length(GND_perm.t_tests(i).clust_info.neg_clust_pval))
        fprintf('     significant: %d\n',length(find(GND_perm.t_tests(i).clust_info.neg_clust_pval<.05)))
        fprintf('     smallest pval: %.3f\n',min(GND_perm.t_tests(i).clust_info.neg_clust_pval))
    end
    
    timepoints = [0:100:1000];
    ntp = length(timepoints);
    xperc = linspace(0.05,0.95,ntp+1);
    
    regidx = [7:12]; %valence correct policy
    
    for iX = regidx
        f = plotTopoRegression(GND_perm,timepoints,iX)
    end
end


%% single-trial regression RESPONSELOCKED

EEGinfix = 'stimLocked\richef-respLocked'; %EEG data for regression
outfix = 'stimLocked\regression'; %folder to save to [outfix model]

%settings for permutation test
timewin = [-650 0]; %-650 is average response time
alpha = 0.05;
nPerm = 10000;

models = {'hierarchicalRL-ar-counter' 'hierarchicalRL-ar' }; 
regtypes = {'resp-revision'}; 
for iM = models
    GNDinfix = strcat(outfix, '-', iM{:});
    for i = 1:length(regtypes)
        EEG_regression(dir,'exp2-eeg',iM{:},regtypes{i},'in',EEGinfix,'out',outfix,'channels',[1:64],'srate',128,'parallel',4); %
        EEG_buildGND(dir,'exp2-eeg',iM{:},GNDinfix,regtypes{i},EEGinfix,[-1 0.5])
        EEG_permuteGND(dir,iM{:},outfix,regtypes{i},EEGinfix,timewin,alpha,nPerm)
    end
end

%plot
if 1
    
    model = 'hierarchicalRL-ar-counter';
    EEGinfix = 'stimLocked\richef-respLocked'; %EEG data for regression
    outfix = 'stimLocked\regression'; %folder to save to [outfix model]
    regtype = 'resp-revision';
    
    targetdir = strcat(outfix, '-', model);
    
    %load results GND
    fin = fullfile(dir.dir_eeg, targetdir, regtype,['GND-' regtype '-perm.GND']);
    load(fin,'-mat');
    
    %plot to get an overview of effects
    for i = 1:length(GND_perm.t_tests)
        fprintf('##########\n')
        fprintf('# %s\n',GND_perm.bin_info(i).bindesc)
        fprintf('##########\n')
        fprintf('regression number: %d\n',i)
        fprintf('   positive clusters: %d\n',length(GND_perm.t_tests(i).clust_info.pos_clust_pval))
        fprintf('     significant: %d\n',length(find(GND_perm.t_tests(i).clust_info.pos_clust_pval<.05)))
        fprintf('     smallest pval: %.3f\n',min(GND_perm.t_tests(i).clust_info.pos_clust_pval))
        fprintf('   negative clusters: %d\n',length(GND_perm.t_tests(i).clust_info.neg_clust_pval))
        fprintf('     significant: %d\n',length(find(GND_perm.t_tests(i).clust_info.neg_clust_pval<.05)))
        fprintf('     smallest pval: %.3f\n',min(GND_perm.t_tests(i).clust_info.neg_clust_pval))
    end
    
    timepoints = [-700:100:200];
    ntp = length(timepoints);
    xperc = linspace(0.05,0.95,ntp+1);
    
    regidx = [2:5]; %valence correct policy
    
    for iX = regidx
        f = plotTopoRegression(GND_perm,timepoints,iX)
    end
    
end



function f = plotTopoRegression(GND_perm,times2plot,regidx)

load('chanlocs.mat','chanlocs')
toposize = 6;

cfg.colors = {'#2D3142','#058ED9','#FF934F','#CC2D35','#E1DAAE', '#848FA2'};
cfg.titlesize = 12;
cfg.infosize = 8;
cfg.alpha = 0.05;

ntp = length(times2plot);
xperc = linspace(0.05,0.95,ntp+1);

for iX = regidx
    
    j = 1;
    f = figure;
    f.Units = 'norm';
%     f.Position = [-1 0.3 0.8 0.2];
    f.Position = [0.1 0.3 0.8 0.2];
    
    titleax = axes('Parent',f);
    titleax.Position = [0.4 0.8 0.2 0.1];
    txt = text(0.5,0.5,sprintf('%s',GND_perm.bin_info(iX).bindesc));
    txt.FontSize = 12;
    txt.FontWeight = 'bold';
    txt.HorizontalAlignment = 'center';
    axis(titleax,'off')
    
    for iT = 1:ntp
        ax(j) = axes('Parent',f);
        ax(j).Position = [xperc(iT) 0.1 1/ntp 0.6];
        
        timept = times2plot(iT); %ms
        timeidx = dsearchn(GND_perm.time_pts',timept);
        dm = squeeze(GND_perm.grands(:,timeidx,iX)); %mean data
        maplimits = [-2 2];
        clustidx = find(timeidx==GND_perm.t_tests(iX).used_tpt_ids);
        chanidx = find(double(GND_perm.t_tests(iX).adj_pval(:,clustidx)<cfg.alpha)==1);
        topoplotMS(dm, chanlocs,'emarker2',{chanidx,'d','k',1},'maplimits',maplimits,'whitebk','on','electrodes','off');
        colormap('redblue')

        txt2 = text(0,0,sprintf('%d ms',times2plot(iT)));
        txt2.FontSize = toposize;
        txt2.HorizontalAlignment = 'center';
        %     txt2.Rotation = 90;
        txt2.Rotation = 0;
        txt2.Position = [0 0.7 0];
        
        gnu = 1;
    end
end





function batch_SurpriseMinimization
% in order to replicate the results from Wurm, Walentowska, Ernst, Severo,
% Pourtois, Steinhauser (2021), download the preprocessed EEG data  and 
% regression results via figshare 
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
dir.setting = 'github'; %adapted for github

%adding model
addpath('Routines\');
addpath('Models\');
addpath('mfit\'); %add Gershman's mfit toolbox

%% Simulation
models = {'hierarchicalRL-full' 'hierarchicalRL-corr' 'hierarchicalRL-incorr'};
SIM_simulateBehavior(models,'empprior','value',100) %model,priors,feedback,simulation

%% Behavioral data
BEH_analysis_regression('exp1-beh') %Fig 2ACE
BEH_analysis_regression('exp2-eeg') %Fig 2BDF

%% Model fitting
optionfit.type = 'behavioral fit';
optionfit.nstarts = [10];
optionfit.rewardscale = [0 100];
optionfit.nBlock = 3;
optionfit.nTrial_learn = 100;
optionfit.nTrial_post = 20;

models = { 'hierarchicalRL-incorr' 'hierarchicalRL-corr'};
for iM = 1:length(models)
    MDL_parafit(dir,'exp1-beh',models{iM},optionfit)
end

models = { 'null' 'jointRL' 'hierarchicalRL-full' 'hierarchicalRL-corr' 'hierarchicalRL-incorr'};
for iM = 1:length(models)
    MDL_parafit(dir,'exp2-eeg',models{iM},optionfit)
end

%% Model comparison
models = {'null' 'jointRL' 'hierarchicalRL-full' 'hierarchicalRL-corr' 'hierarchicalRL-incorr'};
MDL_modelcomp(dir,'exp1-beh',models)
MDL_modelcomp(dir,'exp2-eeg',models)

%% correlate behavior and parameters
MDL_regression(dir,'exp1-beh','hierarchicalRL-full') %Fig 3A
MDL_regression(dir,'exp2-eeg','hierarchicalRL-full') %Fig 3B

%% Simulate variables from fits
optionsim.type = 'simulation';
optionsim.nSamples = 10;
optionsim.fbtype = 'relpay';
optionsim.nSet = 2;
optionsim.nBlock = 3;
optionsim.nTrial_learn = 100; %number of simulated trials
MDL_varsim(dir,'exp2-eeg','hierarchicalRL-full',optionsim)

%% single-trial regression
model = 'hierarchicalRL-full';
EEGinfix = 'fbLocked\richef'; %EEG data for regression
outfix = 'fbLocked\regression'; %folder to save to [outfix model]
GNDinfix = strcat(outfix, '-', model);

% regtypes = {'fb-surprise-stepwise' 'fb-valentin'}; %old: 'hierarchical' 'orthogonal' 'single'
regtypes = {'fb-surprise-stepwise'};
for i = 1:length(regtypes)
    EEG_regression(dir,'exp2-eeg',model,regtypes{i},'in',EEGinfix,'out',outfix,'channels',[1:64],'srate',128,'parallel',0); %
    EEG_buildGND(S,options,GNDinfix,regtypes{i})
end


%% plotting and statistics
model = 'hierarchicalRL-full';
regtype = 'fb-surprise-stepwise';
outfix = 'fbLocked\ERP-regression';
targetdir = strcat(outfix, '-', model);

%load GND (will also be used for ERP plotting)
fin = fullfile(dir.EEGdir, targetdir, regtype,['GND-' regtype '.GND']);
load(fin,'-mat');

%settings for permutation test
timewin = [0 1000]; 
alpha = 0.05;
nPerm = 10000;

% loop over all regressors and save results
GND_perm = GND_rpe;
for i = [1:length(GND_perm.bin_info)]
    fprintf('\n##########\n# %s #\n##########\n',GND_perm.bin_info(i).bindesc)
    GND_perm=clustGND(GND_perm,i,'time_wind',timewin,'thresh_p',.05,'alpha',.05,'n_perm',nPerm,...
        'save_GND','no','plot_raster','no','plot_mn_topo','no');
    gnu = 1;
end
fout = fullfile(dir.EEGdir, targetdir, regtype,['GND-' regtype '-perm.GND']);
% save(fout,'GND_perm');

%load results GND
fin = fullfile(dir.EEGdir, targetdir, regtype,['GND-' regtype '-perm.GND']);
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









function BEH_parafit(dir,fin,model,options)
% S - VP structure
% fn - filename for data
% model - case name for model
% options - specifying environmental and fitting parameters

%% information on model
infopt.type = options.type;
infopt.nstarts = options.nstarts;    % number of random parameter initializations
infopt.rewardscale = options.rewardscale;

%load behavioral data
load(fullfile(dir.dir,['BEHdata-' fin '.mat']),'BEHdata')

%build parameter structure
[~,param,f] = defineModel_eegCAT(model);

VPs = [BEHdata.VPnum];
for i = 1:length(BEHdata)    
    
    r1 = squeeze([BEHdata(i).rew1]);
    r2 = squeeze([BEHdata(i).rew2]);
    r1 = (r1/50)-1;
    r2 = (r2/50)-1;
    
    %build data structure
    data(i).d1 = squeeze([BEHdata(i).dec1])+1;
    data(i).d2 = squeeze([BEHdata(i).dec2])+1;
    data(i).r1 = r1;
    data(i).r2 = r2;
    
    
    data(i).N = options.nTrial_learn*options.nBlock;
    data(i).nB = options.nBlock;
    data(i).nT = options.nTrial_learn;
    data(i).info = infopt;
end

% run optimization
fprintf('... Fitting RL model\n')
fprintf('       Type: %s\n',model);


% results = fw_mfit_optimize2(f,param,data);
% results = mfit_optimize(f,param,data,infopt.nstarts);
results = mfit_optimize_parallel(f,param,data,infopt.nstarts);

fout = fullfile(dir.dir,strcat('Fits\fit-',model,'-',fin));
fprintf('... Saving RL model\n')
fprintf('    File: %s\n',fout)
save(fout,'results')

gnu = 1;

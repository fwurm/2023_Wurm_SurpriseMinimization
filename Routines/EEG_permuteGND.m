function EEG_permuteGND(dir,model,outfix,regtype,EEGinfix,timewin,alpha,nPerm)

%% plotting and statistics
% model = 'hierarchicalRL-counterfact';
% regtype = 'resp-revision'
% outfix = 'stimLocked\regression';
targetdir = strcat(outfix, '-', model);

%load GND (will also be used for ERP plotting)
fin = fullfile(dir.dir_eeg, targetdir, regtype,['GND-' regtype '.GND']);
load(fin,'-mat');


% loop over all regressors and save results
GND_perm = GND_rpe;
for i = [1:length(GND_perm.bin_info)]
    fprintf('\n##########\n# %s #\n##########\n',GND_perm.bin_info(i).bindesc)
    GND_perm=clustGND(GND_perm,i,'time_wind',timewin,'thresh_p',.05,'alpha',.05,'n_perm',nPerm,...
        'save_GND','no','plot_raster','no','plot_mn_topo','no');
    gnu = 1;
end
fout = fullfile(dir.dir_eeg, targetdir, regtype,['GND-' regtype '-perm.GND']);
save(fout,'GND_perm');
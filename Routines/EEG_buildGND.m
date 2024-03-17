function EEG_buildGND(dir,fin,model,inprefix,regtype,eeginfix,framewin)

% settings
% framewin = [-1 0.5]; %time of interest
nChan = 64; %number of channels

%get input folder
filedir = fullfile(dir.dir_eeg, inprefix, regtype);

% get behavior from simluations
[BEH] = EEG_loadVars(dir,fin,model); 
nS = length(BEH);

%prepare regression outputs
for iVP = 1:nS
    
    %load regression data
    fn = fullfile(filedir, sprintf('VP%d.mat',iVP));
    load(fn)
    
    %add to structure
    RPEs(iVP) = data;   
    
    iX = 0; %number of overall coefficients
    for iOut = 1:length(RPEs(1).coefname) %number of regression models
        for iFact = 1:length(RPEs(1).coefname{iOut}) %number of coefficients per model
            iX = iX + 1;
            
            %normalize beta values by their respective SD
            rpevalues = RPEs(iVP).bvals{iOut}(1:nChan,:,iFact);
            rpestd = std(reshape(rpevalues,numel(rpevalues),1));
            RPEdat{iX}(:,:,iVP) = rpevalues./rpestd; %beta normalization
            
            condnames{iX} = ['model' num2str(iOut) '-' strrep(RPEs(1).coefname{iOut}{iFact},':','By')];
            
        end
    end 
       
    %build event struct
    events(iVP).type = 'RPEweight';
    events(iVP).latency = 63.5 + 188*(iVP-1);
    events(iVP).epoch = iVP;
    events(iVP).trials = 'tba';
    events(iVP).setname = ['Regression weight for ' sprintf('VP%d',iVP)];
    
    %build epoch struct
    epochs(iVP).event = iVP;
    epochs(iVP).eventtype = 'RPEweight';
    epochs(iVP).eventlatency = 63.5 + 188*(iVP-1);
    epochs(iVP).eventtrials = 'tba';
    epochs(iVP).eventsetname = ['Regression weight for ' sprintf('VP%d',iVP)];

end
condcount = iX;

%construct timeline
timewin = cell2mat(cellfun(@(x) x(end) - x(1),{RPEs.times},'UniformOutput',false))./1000;
freqrate = timewin./cellfun(@length,{RPEs.times});
srates = floor(1./freqrate);
if length(unique(srates)) > 1
    error('differente sampling rates')
else
    srate = unique(srates);
end
timepnts = 1000.*(framewin(1):1/srate:framewin(2));
times.framewin_sec = framewin;
times.srate = srate;
times.pnts = timepnts;


%load musterGND file (MassUnivariateToolbox)
fn = '.\Routines\musterGND.GND';
load(fn,'-mat');
GND_rpe = musterGND;

%load EEG set (to get chanlocs)
fn = sprintf('VP%d.set',iVP);
raw = pop_loadset('filename',fn,'filepath',fullfile(dir.dir_eeg,eeginfix));
chanlocs = raw.chanlocs(1:nChan);
for i = 1:nChan
    chanlocs(i).X = chanlocs(i).X/100;
    chanlocs(i).Y = chanlocs(i).Y/100;
    chanlocs(i).Z = chanlocs(i).Z/100;
end
raw = pop_select(raw, 'channel', [1:nChan]);
raw_s = pop_resample(raw,srate); %resampling
clear raw

%% Construct EEGlab file 
iX = 0;
for iOut = 1:length(RPEs(1).coefname)
    for iFact = 1:length(RPEs(1).coefname{iOut})
        iX = iX+1;
        
        data = raw_s;
        data.trials = nS;
        data.event = events;
        data.epoch = epochs;
             
        data.setname = ['RegressionWeights: model ' num2str(iOut) ' - ' RPEs(1).coefname{iOut}{iFact}] ;
        data.data = RPEdat{iX};
        
        %save set
        fprintf('Saving RPE regression weights...\n')
        fnout = ['EEG-model' num2str(iOut) '-' strrep(RPEs(1).coefname{iOut}{iFact},':','By') '.set'];
        data = pop_saveset( data, 'filename' , fnout, 'filepath', filedir);
        
    end
end

%% Construct GND for reward prediction errors

fout = ['GND-' regtype];

GND_rpe.exp_desc = 'RPEweights for all VPs';
GND_rpe.filename = fout;
GND_rpe.filepath = filedir;

GND_rpe.grands = [];
GND_rpe.grands_stder = [];
GND_rpe.grands_t = [];
GND_rpe.bin_info = [];
indiv_erps = [];

for i = 1:length(RPEdat)
    grands(:,:,i) = squeeze(mean(RPEdat{i}(:,:,:),3));
    grands_stder(:,:,i) = squeeze(std(RPEdat{i}(:,:,:),[],3)) / sqrt(size(RPEdat{i},3));
    
    GND_rpe.bin_info(i).bindesc = condnames{i};
    GND_rpe.bin_info(i).condcode = 1;
    
    for iVP = 1:nS
        indiv_erps(:,:,i,iVP) = RPEdat{i}(:,:,iVP);
    end
end

grands_t =  grands ./ grands_stder;

GND_rpe.grands = grands;
GND_rpe.grands_stder = grands_stder;
GND_rpe.grands_t = grands_t;
GND_rpe.sub_ct = ones(1,condcount(1))*nS; %VP
GND_rpe.chanlocs = chanlocs;
GND_rpe.time_pts = times.pnts(2:end);
GND_rpe = rmfield(GND_rpe,'bsln_wind');
GND_rpe.srate = srate;
GND_rpe.indiv_fnames = sprintfc('VP%d',1:nS);
GND_rpe.indiv_subnames = sprintfc('VP%d',1:nS);
GND_rpe.indiv_bin_ct = ones(nS,condcount(1));
GND_rpe.indiv_bin_raw_ct = ones(nS,condcount(1));
GND_rpe.indiv_erps = indiv_erps;
GND_rpe.indiv_art_ics = cell(1,nS);

% save GND
fout = [filedir '\GND-' regtype '.GND'];
save(fout,'GND_rpe');
function fw_singletrial_sambrook_gamblearn(dir,fin,model,regtype,varargin)

%define options
options.model = model;
options.regtype = regtype;
options.inprefix = '';
options.outprefix = '';
options.srate = [];
options.channels = [];
options.length_codes = 7; %length of behavioral codes from EEG preprocessing

for iVar = 1:length(varargin)
    if strcmp(varargin{iVar}, 'in')
        if length(varargin)>iVar
            options.inprefix = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''inprefix'' is not valid!');
        end
    end
    if strcmp(varargin{iVar}, 'out')
        if length(varargin)>iVar
            options.outprefix = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''outprefix'' is not valid!');
        end
    end
    if strcmp(varargin{iVar}, 'method')
        if length(varargin)>iVar
            options.method3 = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''beh'' is not valid!');
        end
    end
    if strcmp(varargin{iVar}, 'srate')
        if length(varargin)>iVar
            options.srate = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''srate'' is not valid!');
        end
    end
    if strcmp(varargin{iVar}, 'modelspec')
        if length(varargin)>iVar
            options.modelspec = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''modelspec'' is not valid!');
        end
    end
    if strcmp(varargin{iVar}, 'channels')
        if length(varargin)>iVar
            options.channels = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''channels'' is not valid!');
        end
    end
    if strcmp(varargin{iVar}, 'parallel')
        if length(varargin)>iVar
            options.parallel = varargin{iVar+1};
        else
            disp('ERROR: Input for parameter ''parallel'' is not valid!');
        end
    end
    
end

%load regression models
regmod = EEG_getRegMod(regtype);
options.regmod = regmod;

fprintf('\n Models:\n')
fprintf('\t%s\n',options.regmod.model)

%make folder to save data 
fprintf('Writing output variable...\n')
outpath = fullfile(dir.dir_eeg, [options.outprefix '-' model], regtype);
if ~exist(outpath, 'dir')
    mkdir(outpath)
end
options.outpath = outpath;

% get behavior from simluations
[BEH] = EEG_loadVars(dir,fin,model); 

if options.parallel > 0
    
    %open parallel pool with requested number of workers
    P = gcp('nocreate');
    
    if isempty(P)
        %         npool = feature('numcores');
        P = parpool('local',options.parallel);
    else
        fprintf('Using existing parallel pool with %d workers.\n',P.NumWorkers);
    end

    parfor iS = 1:length(BEH)
        subfun_regress(iS,options,BEH(iS));
    end
else
    for iS = 1:length(BEH)
        subfun_regress(iS,dir,options,BEH(iS));
    end
end

function subfun_regress(iS,dir,options,BEH)

EEGfn = sprintf('VP%d',iS);

t = [];
mdl = [];
bvals = [];
out = struct();

fprintf('#####################\n')
fprintf('# Processing %s... #\n',EEGfn)
fprintf('#####################\n')

%% Read EEG data
inpath = fullfile(dir.dir_eeg, options.inprefix);
raw = pop_loadset('filename',[EEGfn '.set'],'filepath',inpath);

%get trialnumbers, set and sequence/position information
trialnums = [raw.event.originalTrialNumber];
set = [raw.event.type];
set = set(1:options.length_codes:end);
set = str2num(set');
pos = [raw.event.type];
pos = pos(options.length_codes:options.length_codes:end);
pos = str2num(pos');
if (length(pos) ~= length(trialnums)) || (length(set) ~= length(trialnums))
    error('something is off')
end

%apply resampling
if ~isempty(options.srate)
    raw = pop_resample(raw,options.srate);
else
    fprintf('   No resampling applied\n')
end

%get channels
if ~isempty(options.channels)
    raw = pop_select(raw,'channel',options.channels);
else
    fprintf('   No channels deleted\n')
end

%prepare behavioral data for regression
[pre] = EEG_prepVars(BEH,trialnums,set,pos);


%prepare tables for regression
tprep = table();
for iReg = 1:length(options.regmod)
    %add variables from regression model
    outsplit=regexp(options.regmod(iReg).model,'[ +*:()]','split'); %split     
    for iV = 3:length(outsplit)
        if ~isempty(outsplit{iV}) && ~ismember(outsplit{iV},tprep.Properties.VariableNames)
        tprep.(outsplit{iV}) = pre.(outsplit{iV});
        end
    end
    %add variables from splitting conditions
    for iS = 1:length(options.regmod(iReg).identVars)
        if ~ismember(options.regmod(iReg).identVars(iS),tprep.Properties.VariableNames)
        tprep.(options.regmod(iReg).identVars{iS}) = pre.(options.regmod(iReg).identVars{iS});
        end
    end
    %add variables from orthogonalization
    for iO = 1:length(options.regmod(iReg).orthorder)
        if ~ismember(options.regmod(iReg).orthorder(iO),tprep.Properties.VariableNames)
        tprep.(options.regmod(iReg).orthorder{iO}) = pre.(options.regmod(iReg).orthorder{iO});
        end
    end
end

%get EEG
y_raw = double(squeeze([raw.data]));

%create separate orthogonalizes tables
tortho = orthotable(tprep,options.regmod);


%% regression
fprintf('Calculating multiple univariate regressions....\n Channel: ')
for i = 1:raw.nbchan
    fprintf('%s, ',raw.chanlocs(i).labels)
    
    for j = 1:raw.pnts
        
        %get eeg data for channel and timepoint
        y = squeeze(y_raw(i,j,:));
                      
        for k = 1:length(options.regmod)
            
            tnew = tortho{k};
            tnew.('eeg') = y;
            
            nanVars = [options.regmod(k).zVars];
            nanVars{end+1} = 'eeg';
            tnew = nantable(tnew,nanVars); %remove nan from table
            
            tnew = identtable(tnew,options.regmod(k).identVars,options.regmod(k).identVals,0); % index table (find specific condition as specified in regmod
                
            tnew = ztable(tnew,options.regmod(k).zVars); %zscore relevant columns           
            
            if strcmp(options.regtype,'choice') 
                tnew = rttable(tnew,rts); %apply reaction time trimming
            end
            
            mdl{k} = fitglm(tnew,options.regmod(k).model,'Distribution','normal'); %fit regression
            bvals{k}(i,j,:) = table2array(mdl{k}.Coefficients(:,1)); %save beta values
        end        
        gnu = 1;
    end
    gnu = 1;
end
fprintf('\n')

%save names of coefficients
coefname = [];
for m = 1:length(options.regmod)
    coefname{m} = mdl{m}.CoefficientNames;
    coefname{m} = cellfun(@(x) regexprep(x,'(\(|\))',''),coefname{m},'UniformOutput',false);
end

%% construct and write output
fprintf('Constructing variable with regression values...\n')

out.VPname = sprintf('VP%03.f',iS);
out.VPnum = iS;
out.regmod = options.regmod;
out.model = options.model;
out.bvals = bvals;
out.coefname = coefname;
out.times = raw.times;
out.options = options;

fprintf('Writing output variable...\n')
fout = fullfile(options.outpath, sprintf('VP%d.mat',iS));
parsave(fout,out);

gnu = 1;


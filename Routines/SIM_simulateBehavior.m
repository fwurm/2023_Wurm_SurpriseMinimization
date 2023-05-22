function simulateBehavior(models,paratype,fbtype,nSubjects)
% input parameter
%   - models
%   - paratype (
%       empprior = parameters drawn from empirical prior distributions
%       random = drawn from random (=uniform) distributions
%       fitted = drawn from behavioral data
%   - fbtype (
%       value = magnitude (reward ranging between 0 and 100)
%       probs = reward probability (reward is either 0 or 1)


%% setting up the environment (identical to experiment)
para_sim.nB = 3; % number of blocks
para_sim.nS = 2; % number of sets per block
para_sim.nT = 100; % number of trials
% rng(0,'twister'); %initialize the random number generator for replicability
para_sim.lowerBound = 1;
para_sim.upperBound = 99;
para_sim.walkMean = 0;
para_sim.walkStd = 20;

walk = struct();
walk.info.type = 'simulation';
walk.info.fbtype = fbtype;
walk.N = para_sim.nT*para_sim.nB*2;
walk.nB = para_sim.nB;
walk.nT = para_sim.nT;


%% retrieve parameters for simulation
if ismember(paratype,{'empprior','random'})   
    for iModel = 1:length(models)
        [mstruct(iModel), ~, ~] = defineModel_eegCAT(models{iModel});
        paraspace{iModel} = sampleParameters(mstruct(iModel),nSubjects,paratype);
    end   
elseif strcmp(paratype,'fitted')
    for iModel = 1:length(models)
        fn_estimation = [dir.dir_fits 'fit_subject_' models{iModel} '_all_v4.mat'];
        load(fn_estimation)
        behavior.nVP = size(results.x,1);
        paraspace{iModel} = [results.x];
    end
else
    error('no suffix for parameter recovery specified');
end

% warning('Parameterspace has been set to equal')
% paraspace{2} = paraspace{1}(:,1:3);
% paraspace{3} = paraspace{1}(:,1:3);

%% generating behavior
fprintf('Starting simulation...\n')
fprintf('Agent 1.')
for i = 1:nSubjects
    
    if rem(i,100)==0
        fprintf('%d.',i)
    end
    
    % generate walks
    [rew_payoff,rew_dichot] = doWalk(para_sim);
    if strcmp(walk.info.fbtype,'value')
        walk.r1(:,:,1) = squeeze(rew_payoff(:,1,:));
        walk.r1(:,:,2) = 100 - squeeze(rew_payoff(:,1,:));
        walk.r2(:,:,1) = squeeze(rew_payoff(:,2,:));
        walk.r2(:,:,2) = 100 - squeeze(rew_payoff(:,2,:));
        walk.r1 = (walk.r1./50)-1;
        walk.r2 = (walk.r2./50)-1;
    elseif strcmp(walk.info.fbtype,'probs')
        walk.r1(:,:,1) = squeeze(rew_dichot(:,1,:));
        walk.r1(:,:,2) = 1 - squeeze(rew_dichot(:,1,:));
        walk.r2(:,:,1) = squeeze(rew_dichot(:,2,:));
        walk.r2(:,:,2) = 1 - squeeze(rew_dichot(:,2,:));
        walk.r1 = (walk.r1.*2)-1;
        walk.r2 = (walk.r2.*2)-1;
    end
    
    %run models
    for iModel = 1:length(models)
        parasubspace = paraspace{iModel}(:,i);
        
        if strcmp(models{iModel},'hierarchicalRL-corr')
            parasubspace(4) = 999; %correct mapping
        elseif strcmp(models{iModel},'hierarchicalRL-incorr')
            parasubspace(4) = -999; %incorrect mapping
        end
        
        data_sim = mstruct(iModel).funhandle(parasubspace,walk);
        
        BEHdata{iModel}(i).VPname = ['Sim' num2str(i)];
        BEHdata{iModel}(i).VPnum = i;
        BEHdata{iModel}(i).dec1 = data_sim.c1;
        BEHdata{iModel}(i).rew1 = data_sim.r1;
        BEHdata{iModel}(i).dec2 = data_sim.c2;
        BEHdata{iModel}(i).rew2 = data_sim.r2;
        
        performance(i,iModel) = mean([mean(mean(data_sim.corr1,2),1) mean(mean(data_sim.corr2,2),1)]);
        
        if strcmp(models{iModel},'hierarchicalRL-corr')
            %                 [pd] = fitdist(reshape(data_sim.RPEp1,1,numel(data_sim.RPEp1))','Normal');
            surprise(i,iModel) = mean(abs(reshape(data_sim.RPEp1,1,numel(data_sim.RPEp1))));
        elseif strcmp(models{iModel},'hierarchicalRL-incorr')
            %                 [pd] = fitdist(reshape(data_sim.RPEp2,1,numel(data_sim.RPEp2))','Normal');
            surprise(i,iModel) = mean(abs(reshape(data_sim.RPEp2,1,numel(data_sim.RPEp2))));
        elseif strcmp(models{iModel},'hierarchicalRL-full')
            cat(:,:,i) = [data_sim.cat];
            %                 warning()
        end
        
        gnu = 1;
    end
    
    gnu = 1;
    
    
    
end
fprintf('\n')



model = 'hierarchicalRL-corr';
fig3A = calcswitch(model,BEHdata{ismember(models,model)}(1:30));

model = 'hierarchicalRL-incorr';
fig3B = calcswitch(model,BEHdata{ismember(models,model)}(1:30));

model = 'hierarchicalRL-full';
fig3C = calcswitch(model,BEHdata{ismember(models,model)}(1:30));

Match=cellfun(@(x) ismember(x, {'hierarchicalRL-corr' 'hierarchicalRL-incorr'}), models, 'UniformOutput', 0);
mapident = cell2mat(Match);

% figure 3D
fig3D = figure;
ax3D = axes;
hist(surprise(:,mapident),100)
xlabel('Surprise')
ylabel('Frequency')
legend({'Correct policy' 'Incorrect policy'})
set(gcf,'color','w');
title('Surprise distribution')

% figure 3E
fig3E = figure;
ax3E = axes;
hold on
plot([1:100],cat(1,:,1),'LineWidth',2)
plot([1:100],cat(2,:,1),'LineWidth',2)
plot([1:100],cat(3,:,1),'LineWidth',2)
plot([1:100],mean(cat,[1 3]),'k','LineWidth',2)
legend({'Block 1' 'Block 2' 'Block 3' 'Mean'},'Location','southeast')
title('Evidence accumulation')




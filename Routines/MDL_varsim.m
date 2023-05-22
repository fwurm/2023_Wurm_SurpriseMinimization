function fw_simBEH(dir,fin,model,options)

%load behavioral data
load(fullfile(dir.dir,strcat('BEHdata-',fin,'.mat')),'BEHdata');

%load fits
load(fullfile(dir.dir_fits,strcat('fit-', model, '-', fin, '.mat')),'results');

%build parameter structure (taken from Gershman example)
[~,~,f] = defineModel_eegCAT(model);


%simulate data based on fits
for i = 1:length(BEHdata)
    
    %variables to simulate
    rpe_p1 = [];
    Qupd_p1 = [];
    rpe_p2 = [];
    Qupd_p2 = [];
    HPE = [];
    
    %random walk info
    walk = [];
    walk.d1 = [BEHdata(i).dec1 + 1];
    walk.d2 = [BEHdata(i).dec2 + 1];
    walk.info = options;
    walk.N = options.nTrial_learn;
    walk.nB = options.nBlock;
    walk.nT = options.nTrial_learn;
    
    %rescale feedback
    if strcmp(options.fbtype,'relpay') %center around 0, with [-1 1]
        walk.r1 = ([BEHdata(i).rew1]./50)-1;
        walk.r2 = ([BEHdata(i).rew2]./50)-1;
    elseif strcmp(options.fbtype,'abs') %center around 0, with [-1 1]
        walk.r1 = ([BEHdata(i).rew1].*2)-1;
        walk.r2 = ([BEHdata(i).rew2].*2)-1;
    else
        warning('no rewardscale specified')
    end

    %get fit parameter
    param = squeeze([results.x(i,:)]);
    
    %run model with parameters
    data = f(param,walk);
        
    c1 = reshape(permute([data.c1],[2,1]),[300 1]);
    c2 = reshape(permute([data.c2],[2,1]),[300 1]);
    
    %policy 1 (correct)
    rpe_p1 = reshape(permute([data.RPEp1],[2,1,3]),[300,2]);
    Q1_p1 = reshape(permute([data.v1p1_upd],[2,1,3]),[300,2]); %[data.v1p1_upd(:,:,[data.c1])];
    Q2_p1 = reshape(permute([data.v2p1_upd],[2,1,3]),[300,2]); %[data.v2p1_upd(:,:,[data.c2])];
%     Qupd_p1 = reshape(permute(Qupd_p1,[2,1,3]),[300,2]);
    
    %policy 1 (incorrect)
    rpe_p2 = reshape(permute([data.RPEp2],[2,1,3]),[300,2]);
    Q1_p2 = reshape(permute([data.v1p2_upd],[2,1,3]),[300,2]); %[data.v1p2_upd(:,:,[data.c1])];
    Q2_p2 = reshape(permute([data.v2p2_upd],[2,1,3]),[300,2]); %[data.v2p2_upd(:,:,[data.c2])];
%     Qupd_p2 = reshape(permute(Qupd_p2,[2,1,3]),[300,2]);
    
    %hierarchical policy
    HPE = reshape([data.HPE]',[300,1]);
    rpe_comp = reshape(permute([data.rpe_comp],[2,1,3]),[300,2]);
    cat = reshape(permute([data.cat],[2,1]),[300 1]);
    cat_upd = reshape(permute([data.cat_upd],[2,1]),[300 1]);
        
    %write variables
    RPEs(i).VPname = sprintf('sub%03.f',i);
    RPEs(i).VPnum = i;
    RPEs(i).c1 = c1;
    RPEs(i).c2 = c2;
    RPEs(i).RPEs_p1 = rpe_p1;
    RPEs(i).Q1_p1 = Q1_p1;
    RPEs(i).Q2_p1 = Q2_p1;
    RPEs(i).RPEs_p2 = rpe_p2;
    RPEs(i).Q1_p2 = Q1_p2;
    RPEs(i).Q2_p2 = Q2_p2;
    RPEs(i).HPE = HPE;
    RPEs(i).RPE_comp = rpe_comp;
    RPEs(i).cat = cat;
    RPEs(i).cat_upd = cat_upd;
    RPEs(i).param = param;
            
    gnu = 1;
end

%save variables
fout = fullfile(dir.dir_fits,strcat('sim-', model, '-', fin, '.mat'));
save(fout,'RPEs')
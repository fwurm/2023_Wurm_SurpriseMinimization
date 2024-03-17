function [BEH] = getBehavior_eegCAT(dir,fin,model)

%% Read Behavioral data
filename = fullfile(dir.dir, strcat('BEHdata-',fin, '.mat'));
fprintf('   Load File: %s\n', filename);
load(filename,'BEHdata');

filename2 = fullfile(dir.dir_fits ,strcat('sim-', model, '-',fin, '.mat'));
fprintf('   Load File: %s\n', filename2);
load(filename2,'RPEs');

value_Q1p1 = cellfun(@(x) diff(x,[],2),{RPEs.Q1_p1},'UniformOutput',false);
value_Q1p2 = cellfun(@(x) diff(x,[],2),{RPEs.Q1_p2},'UniformOutput',false);
diff_Q1p1 = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q1_p1},'UniformOutput',false);
diff_Q1p2 = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q1_p2},'UniformOutput',false);

value_Q2p1 = cellfun(@(x) diff(x,[],2),{RPEs.Q2_p1},'UniformOutput',false);
value_Q2p2 = cellfun(@(x) diff(x,[],2),{RPEs.Q2_p2},'UniformOutput',false);
diff_Q2p1 = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q2_p1},'UniformOutput',false);
diff_Q2p2 = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q2_p2},'UniformOutput',false);

value_Q1p1_old = cellfun(@(x) diff(x,[],2),{RPEs.Q1_p1_old},'UniformOutput',false);
value_Q1p2_old = cellfun(@(x) diff(x,[],2),{RPEs.Q1_p2_old},'UniformOutput',false);
diff_Q1p1_old = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q1_p1_old},'UniformOutput',false);
diff_Q1p2_old = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q1_p2_old},'UniformOutput',false);

value_Q2p1_old = cellfun(@(x) diff(x,[],2),{RPEs.Q2_p1_old},'UniformOutput',false);
value_Q2p2_old = cellfun(@(x) diff(x,[],2),{RPEs.Q2_p2_old},'UniformOutput',false);
diff_Q2p1_old = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q2_p1_old},'UniformOutput',false);
diff_Q2p2_old = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q2_p2_old},'UniformOutput',false);

value_Q1net_old = cellfun(@(x) diff(x,[],2),{RPEs.Q1net_old},'UniformOutput',false);
value_Q2net_old = cellfun(@(x) diff(x,[],2),{RPEs.Q2net_old},'UniformOutput',false);
diff_Q1net_old = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q1net_old},'UniformOutput',false);
diff_Q2net_old = cellfun(@(x) abs(diff(x,[],2)),{RPEs.Q2net_old},'UniformOutput',false);

BEH = [];
for i = 1:length([RPEs.VPnum])
    BEH(i).RPEvpnums = [RPEs(i).VPnum];
    BEH(i).BEHvpnums = [BEHdata(i).VPnum];
    
    
    BEH(i).c1 = [RPEs(i).c1];
    BEH(i).c2 = [RPEs(i).c2];
    
    BEH(i).rpe_p1 = [RPEs(i).RPEs_p1];
    BEH(i).rpe_p2 = [RPEs(i).RPEs_p2];
    
    BEH(i).Q1_p1 = [RPEs(i).Q1_p1];
    BEH(i).Q1_p2 = [RPEs(i).Q1_p2];
    BEH(i).value1_p1 = value_Q1p1{i};
    BEH(i).value1_p2 = value_Q1p2{i};
    BEH(i).diff1_p1 = diff_Q1p1{i};
    BEH(i).diff1_p2 = diff_Q1p2{i};
    
    BEH(i).Q2_p1 = [RPEs(i).Q2_p1];
    BEH(i).Q2_p2 = [RPEs(i).Q2_p2];
    BEH(i).value2_p1 = value_Q2p1{i};
    BEH(i).value2_p2 = value_Q2p2{i};
    BEH(i).diff2_p1 = diff_Q2p1{i};
    BEH(i).diff2_p2 = diff_Q2p2{i};
    
    BEH(i).Q1_p1_old = [RPEs(i).Q1_p1_old];
    BEH(i).Q1_p2_old = [RPEs(i).Q1_p2_old];
    BEH(i).value1_p1_old = value_Q1p1_old{i};
    BEH(i).value1_p2_old = value_Q1p2_old{i};
    BEH(i).diff1_p1_old = diff_Q1p1_old{i};
    BEH(i).diff1_p2_old = diff_Q1p2_old{i};
    
    BEH(i).Q2_p1_old = [RPEs(i).Q2_p1_old];
    BEH(i).Q2_p2_old = [RPEs(i).Q2_p2_old];
    BEH(i).value2_p1_old = value_Q2p1_old{i};
    BEH(i).value2_p2_old = value_Q2p2_old{i};
    BEH(i).diff2_p1_old = diff_Q2p1_old{i};
    BEH(i).diff2_p2_old = diff_Q2p2_old{i};
    
    BEH(i).value_Q1net_old = value_Q1net_old{i};
    BEH(i).value_Q2net_old = value_Q2net_old{i};
    BEH(i).diff_Q1net_old = diff_Q1net_old{i};
    BEH(i).diff_Q2net_old = diff_Q2net_old{i};
    
    BEH(i).evidence_fb = [RPEs(i).evidence_fb]; %evidence
    BEH(i).evidence_resp = [RPEs(i).evidence_resp]; %evidence
    BEH(i).caw_old = [RPEs(i).cat]; %credit assignmet weight
    BEH(i).caw_new = [RPEs(i).cat_upd]; %updated credit assignmet weight
    
    BEH(i).param = RPEs(i).param;
    
end

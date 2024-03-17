function regmod = getRegMod(regtype)


i = 0; %number of models

if strcmp(regtype,'fb-surprise-stepwise')
    i = i+1; %standard RL model
    regmod(i).name = 'splitRPE-only';
    regmod(i).model = 'eeg ~ (valenceP1 * surpriseP1) + (valenceP2 * surpriseP2)';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'valenceP1' 'surpriseP1' 'valenceP2' 'surpriseP2'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P1-positive';
    regmod(i).model = 'eeg ~ surpriseP1';
    regmod(i).identVars = {'valenceP1'}; regmod(i).identVals = {1};
    regmod(i).zVars = {'surpriseP1'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P1-negative';
    regmod(i).model = 'eeg ~ surpriseP1';
    regmod(i).identVars = {'valenceP1'}; regmod(i).identVals = {-1};
    regmod(i).zVars = {'surpriseP1'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P2-positive';
    regmod(i).model = 'eeg ~ surpriseP2';
    regmod(i).identVars = {'valenceP2'}; regmod(i).identVals = {1};
    regmod(i).zVars = {'surpriseP2'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P2-negative';
    regmod(i).model = 'eeg ~ surpriseP2';
    regmod(i).identVars = {'valenceP2'}; regmod(i).identVals = {-1};
    regmod(i).zVars = {'surpriseP2'};
    regmod(i).orthorder = {};
    
    i = i+1; %hierarchical evidence model
    regmod(i).name = 'HPE-only';
    regmod(i).model = 'eeg ~ rpeP1 + rpeP2 + hpe';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'rpeP1' 'rpeP2' 'hpe'};
    regmod(i).orthorder = {};
    
elseif strcmp(regtype,'fb-revision')
    
    i = i+1; %standard RL model
    regmod(i).name = 'splitRPE-only';
    regmod(i).model = 'eeg ~ rpeP1 + rpeP2 + evidence_fb + caw';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'rpeP1' 'rpeP2' 'evidence_fb' 'caw'};
    regmod(i).orthorder = {};
    
    i = i+1; %standard RL model
    regmod(i).name = 'splitRPE-only';
    regmod(i).model = 'eeg ~ (valenceP1 * surpriseP1) + (valenceP2 * surpriseP2)';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'valenceP1' 'surpriseP1' 'valenceP2' 'surpriseP2'};
    regmod(i).orthorder = {};   
    
elseif strcmp(regtype,'fb-revision-posthoc')
    
    i = i+1;
    regmod(i).name = 'P1-positive';
    regmod(i).model = 'eeg ~ surpriseP1';
    regmod(i).identVars = {'valenceP1'}; regmod(i).identVals = {1};
    regmod(i).zVars = {'surpriseP1'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P1-negative';
    regmod(i).model = 'eeg ~ surpriseP1';
    regmod(i).identVars = {'valenceP1'}; regmod(i).identVals = {-1};
    regmod(i).zVars = {'surpriseP1'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P2-positive';
    regmod(i).model = 'eeg ~ surpriseP2';
    regmod(i).identVars = {'valenceP2'}; regmod(i).identVals = {1};
    regmod(i).zVars = {'surpriseP2'};
    regmod(i).orthorder = {};
    
    i = i+1;
    regmod(i).name = 'P2-negative';
    regmod(i).model = 'eeg ~ surpriseP2';
    regmod(i).identVars = {'valenceP2'}; regmod(i).identVals = {-1};
    regmod(i).zVars = {'surpriseP2'};
    regmod(i).orthorder = {};
    

elseif strcmp(regtype,'fb-startend')
    
    i = i+1; %standard RL model
    regmod(i).name = 'splitRPE-startend';
    regmod(i).model = 'eeg ~ startend * ((valenceP1 * surpriseP1) + (valenceP2 * surpriseP2) + surpriseP1:surpriseP2)';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'valenceP1' 'surpriseP1' 'valenceP2' 'surpriseP2' 'startend'};
    regmod(i).orthorder = {};
    
    i = i+1; %standard RL model
    regmod(i).name = 'RPE-startend';
    regmod(i).model = 'eeg ~ startend * (rpeP1 + rpeP2 + evidence_fb + caw)';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'rpeP1' 'rpeP2' 'evidence_fb' 'caw'};
    regmod(i).orthorder = {};    
    
    
elseif strcmp(regtype,'resp-revision')
    
    i = i+1; %standard RL model
    regmod(i).name = 'GLM4-evidence1';
    regmod(i).model = 'eeg ~ rpeP1_old + rpeP2_old + hpe_old + caw_old ';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'rpeP1_old' 'rpeP2_old' 'hpe_old' 'caw_old'};
    regmod(i).orthorder = {};

else
    error('unknown regression model')
end


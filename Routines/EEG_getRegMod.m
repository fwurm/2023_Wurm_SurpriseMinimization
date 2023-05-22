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
    regmod(i).model = 'eeg ~ hpe';
    regmod(i).identVars = {}; regmod(i).identVals = {};
    regmod(i).zVars = {'hpe'};
    regmod(i).orthorder = {};
    
else
    error('unknown regression model')
end


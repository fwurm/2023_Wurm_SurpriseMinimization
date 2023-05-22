function [mstruct, param,f] = defineModel(modelname)
%build structure mstruct, which contains the hard-coded parameters in this
%function according to string input modelname


if strcmp(modelname,'jointRL')
    mstruct.modelname = 'Joint-action RL model with perseveration';
    mstruct.modelshort = 'jointRL';
    mstruct.funhandle = @model_eegCAT_jointRL;
    
    mstruct.name = {'learning rate' 'inverse temperature' 'perseveration'};          
    a1 = 1.2; a2 = 1.2; %learning rate   
    b1 = 2;  b2 = 1; %inverse temperature
    e1 = 0; e2 = 1; %perseveration
    
    lb = [0 0 -5];
    ub = [1 20 5];
    samplehandle_empprior = @(x) [betarnd(a1,a2,x,1) gamrnd(b1,b2,x,1) normrnd(d1,d2,x,1) ];        

    param(1).name = 'learning rate';
    param(1).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
    param(1).lb = 0; 
    param(1).ub = 1;

    param(2).name = 'inverse temperature';
    param(2).logpdf = @(x) sum(log(gampdf(x,b1,b2)));  % log density function for prior
    param(2).lb = 0; 
    param(2).ub = 20; 

    param(3).name = 'perseveration';
    param(3).logpdf = @(x) sum(log(normpdf(x,e1,e2)));
    param(3).lb = -5; 
    param(3).ub = 5;
    
    f = @model_eegCAT_jointRL;  
    
elseif strcmp(modelname,'hierarchicalRL-full')
    mstruct.modelname = 'RL model with surprise minimization and perseveration';
    mstruct.modelshort = 'hierarchicalRL';
    mstruct.funhandle = @model_eegCAT_hierarchicalRL;
    
    mstruct.name = {'learning rate' 'inverse temperature' 'perseveration' 'assignment rate'};          
    a1 = 1.2; a2 = 1.2; %learning rate   
    b1 = 2;  b2 = 1; %inverse temperature
    e1 = 0; e2 = 1; %perseveration
    a1 = 1.2; a2 = 1.2; %assignment rate  
    
    lb = [0 0 -5 0];
    ub = [1 20 5 1];
    samplehandle_empprior = @(x) [betarnd(a1,a2,x,1) gamrnd(b1,b2,x,1) normrnd(e1,e2,x,1) betarnd(a1,a2,x,1)];        

    param(1).name = 'learning rate';
    param(1).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
    param(1).lb = 0; 
    param(1).ub = 1;

    param(2).name = 'inverse temperature';
    param(2).logpdf = @(x) sum(log(gampdf(x,b1,b2)));  % log density function for prior
    param(2).lb = 0; 
    param(2).ub = 20; 

    param(3).name = 'perseveration';
    param(3).logpdf = @(x) sum(log(normpdf(x,e1,e2)));
    param(3).lb = -5; 
    param(3).ub = 5;
    
    param(4).name = 'assignment rate';
    param(4).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
    param(4).lb = 0; 
    param(4).ub = 1;
    
    f = @model_eegCAT_hierarchicalRL; 

elseif strcmp(modelname,'hierarchicalRL-corr') || strcmp(modelname,'hierarchicalRL-incorr')
    mstruct.modelname = 'RL model with surprise minimization and perseveration and fixed mapping';
    mstruct.modelshort = modelname;
    mstruct.funhandle = @model_eegCAT_hierarchicalRL;
    
    mstruct.name = {'learning rate' 'inverse temperature' 'perseveration' 'assignment rate'};          
    a1 = 1.2; a2 = 1.2; %learning rate   
    b1 = 2;  b2 = 1; %inverse temperature
    e1 = 0; e2 = 1; %perseveration
    
    lb = [0 0 -5];
    ub = [1 20 5];
    samplehandle_empprior = @(x) [betarnd(a1,a2,x,1) gamrnd(b1,b2,x,1) normrnd(e1,e2,x,1)];        

    param(1).name = 'learning rate';
    param(1).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
    param(1).lb = 0; 
    param(1).ub = 1;

    param(2).name = 'inverse temperature';
    param(2).logpdf = @(x) sum(log(gampdf(x,b1,b2)));  % log density function for prior
    param(2).lb = 0; 
    param(2).ub = 20; 

    param(3).name = 'perseveration';
    param(3).logpdf = @(x) sum(log(normpdf(x,e1,e2)));
    param(3).lb = -5; 
    param(3).ub = 5;
    
    if strcmp(modelname,'hierarchicalRL-corr')
        f = @model_eegCAT_hierarchicalRL_corr;
    elseif strcmp(modelname,'hierarchicalRL-incorr')
        f = @model_eegCAT_hierarchicalRL_incorr;
    end
    
elseif strcmp(modelname,'null')
    mstruct.modelname = 'Null model';
    mstruct.modelshort = 'null';
    mstruct.funhandle = @model_eegCAT_null;
    
    mstruct.name = {'inverse temperature'};          
    b1 = 2;  b2 = 1; %inverse temperature
    
    lb = [0];
    ub = [20];
    samplehandle_empprior = @(x) [gamrnd(b1,b2,x,1)];        

    param(1).name = 'inverse temperature';
    param(1).logpdf = @(x) sum(log(gampdf(x,b1,b2)));  % log density function for prior
    param(1).lb = 0; 
    param(1).ub = 20; 
    
    f = @model_eegCAT_null;     
   
    
elseif strcmp(modelname,'hypothesisRL')
    mstruct.modelname = 'RL model with hypothesis testing mechanisms and perseveration';
    mstruct.modelshort = 'hypothesisRL';
    mstruct.funhandle = @model_eegCAT_hypothesisRL;
    
    mstruct.name = {'learning rate' 'inverse temperature' 'perseveration' 'assignment rate'};          
    a1 = 1.2; a2 = 1.2; %learning rate   
    b1 = 2;  b2 = 1; %inverse temperature
    e1 = 0; e2 = 1; %perseveration
    f1 = exp(1); f2 = exp(1); %assignment rate  
    
%     lb = [0 0 0.5 0 -1];
%     ub = [1 5 1 0.5 1];
    lb = [0 0 -5 0 0];
    ub = [1 20 5 1 5];
    samplehandle_empprior = @(x) [betarnd(a1,a2,x,1) gamrnd(b1,b2,x,1) normrnd(d1,d2,x,1) betarnd(a1,a2,x,1) gamrnd(f1,f2,x,1) ];        

    param(1).name = 'learning rate';
    param(1).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
    param(1).lb = 0; 
    param(1).ub = 1;

    param(2).name = 'inverse temperature';
    param(2).logpdf = @(x) sum(log(gampdf(x,b1,b2)));  % log density function for prior
    param(2).lb = 0; 
    param(2).ub = 20; 

    param(3).name = 'perseveration';
    param(3).logpdf = @(x) sum(log(normpdf(x,e1,e2)));
    param(3).lb = -5; 
    param(3).ub = 5;
    
    a = 1.2; b = 1.2;  % parameters of the gamma prior
    param(4).name = 'evidence decay';
    param(4).logpdf = @(x) sum(log(betapdf(x,a1,b1)));  % log density function for prior
    param(4).lb = 0;    % lower bound
    param(4).ub = 1;   % upper bound
    
    a = exp(1); b = exp(1);  % parameters of the gamma prior
    param(5).name = 'decision threshold';
    param(5).logpdf = @(x) sum(log(gampdf(x,f1,f2)));  % log density function for prior
    param(5).lb = 0;    % lower bound
    param(5).ub = 5;   % upper bound 
    
    f = @model_eegCAT_hypothesisRL; 
    
else
    warning('unknown model specification - model not coded')    
end

mstruct.lb = lb;
mstruct.ub = ub;
mstruct.samplehandle_empprior = samplehandle_empprior;
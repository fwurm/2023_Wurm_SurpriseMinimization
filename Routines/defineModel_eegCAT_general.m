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
    samplehandle_empprior_name = '@(x) [betarnd(a1,a2,x,1) gamrnd(b1,b2,x,1) normrnd(e1,e2,x,1)';        

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
    
elseif contains(modelname,'hierarchicalRL')
    mstruct.modelname = 'RL model with surprise minimization and perseveration';
    mstruct.modelshort = 'hierarchicalRL';
    mstruct.funhandle = @model_eegCAT_hierarchicalRL_general;
    
    mstruct.name = {'learning rate' 'inverse temperature' 'perseveration'};          
    a1 = 1.2; a2 = 1.2; %learning rate   
    b1 = 2;  b2 = 1; %inverse temperature
    e1 = 0; e2 = 1; %perseveration
    a1 = 1.2; a2 = 1.2; %assignment rate  
    e1 = 0; e2 = 1; %starting point
    a1 = 1.2; a2 = 1.2; %learning decay
    
    lb = [0 0 -5];
    ub = [1 20 5];
    samplehandle_empprior_name = '@(x) [betarnd(a1,a2,x,1) gamrnd(b1,b2,x,1) normrnd(e1,e2,x,1)';        

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
    
    iPara = 3;
    if contains(modelname,'-ar')
        iPara = iPara + 1; 
        mstruct.name{iPara} = 'assignment rate';
        lb = [lb 0];
        ub = [ub 1];
        samplehandle_empprior_name = cat(2,samplehandle_empprior_name,' betarnd(a1,a2,x,1)');
        
        param(iPara).name = 'assignment rate';
        param(iPara).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
        param(iPara).lb = 0; 
        param(iPara).ub = 1;       
    end
    
    if contains(modelname,'-start')
        iPara = iPara + 1; 
        mstruct.name{iPara} = 'start point';
        lb = [lb -5];
        ub = [ub 5];
        samplehandle_empprior_name = cat(2,samplehandle_empprior_name,' normrnd(e1,e2,x,1)');
        
        param(iPara).name = 'starting point';
        param(iPara).logpdf = @(x) sum(log(normpdf(x,e1,e2)));
        param(iPara).lb = -5; 
        param(iPara).ub = 5;      
    end
    
    if contains(modelname,'-decay')
        iPara = iPara + 1; 
        mstruct.name{iPara} = 'learning decay';
        lb = [lb 0];
        ub = [ub 1];
        samplehandle_empprior_name = cat(2,samplehandle_empprior_name,' betarnd(a1,a2,x,1)');
        
        param(iPara).name = 'learning decay';
        param(iPara).logpdf = @(x) sum(log(betapdf(x,a1,a2)));
        param(iPara).lb = 0; 
        param(iPara).ub = 1;       
    end
    
    
       
    f = @model_eegCAT_hierarchicalRL_general; 
    

     
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
    
else
    warning('unknown model specification - model not coded')    
end

mstruct.lb = lb;
mstruct.ub = ub;
mstruct.samplehandle_empprior = eval(cat(2,samplehandle_empprior_name, ']'));

gnu = 1;
function out = sampleParameters(mstruct,nSamples,option)

out = [];


if strcmp(option,'empprior')
    display('sampling from empirical priors...')
    
    for i = 1:nSamples       
        draw = mstruct.samplehandle_empprior(1);
        while any(draw>mstruct.ub) || any(draw<mstruct.lb)
            draw = mstruct.samplehandle_empprior(1);
        end
        out(:,i) = draw;
    end
    
    
elseif strcmp(option,'random')
    display('sampling from uniform priors...')
     
    for iParam = 1:length(mstruct.lb)
        out(iParam,:) = unifrnd(mstruct.lb(iParam),mstruct.ub(iParam),nSamples,1);
    end
    
else
    error('check input ''option''')
end
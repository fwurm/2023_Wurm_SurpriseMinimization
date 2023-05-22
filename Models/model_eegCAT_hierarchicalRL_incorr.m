function output = model_postCAT(x,data)

% Likelihood function for Q-learning on experiment eegCATrandom with
% parameters (1)learning rate, (2)inverse temperature, (3)
% perseveration (4)assignment rate
%
% USAGE: lik = model_postCAT(x,data)
%
% INPUTS:
%   x - parameters:
%       x(1) - learning rate
%       x(2) - inverse temperature
%       x(3) - perseveration
%       x(4) - assignment rate
%   data - structure with the following fields
%          .d1 - [1 x nT] choices/decisions stage1
%          .d2 - [1 x nT] choices/decisions stage2
%          .r(1) - [1 x nT] rewards stage1
%          .r(2) - [1 x nT] rewards stage2
%          .N - number of trials for experiment
%          .nB - number of blocks
%          .nT - number of trials
%
% OUTPUTS:
%   lik - log-likelihood
%
% Franz Wurm, October 2021
% adapted from Sam Gershman, June 2015

%get parameters
lr = x(1);
it = x(2);
pe = x(3);
ar = 0;
sb = -100;

if isfield(data,'d1') && isfield(data,'d2')
    actiontype = 'decision_given';
else
    actiontype = 'sample_decision';
end

info = [data.info];

C = [2 2]; %two options for each action

%initialize
lik = 0; %likelihood (sum over log)

for iB = 1:data.nB
    
    v1p1 = zeros(1,C(1));  % initial values stage 1 policy 1
    v2p1 = zeros(1,C(2));  % initial values stage 2 policy 1
    v1p2 = zeros(1,C(1));  % initial values stage 1 policy 2
    v2p2 = zeros(1,C(2));  % initial values stage 2 policy 2
    
    ca = sb; % initialize credit assignment
    cat = exp(ca)/(1+exp(ca)); %rescale to [0 1] using inverse logit
    
    %initialize indicator for perseveration
    ind1 = zeros(1,C(1)); 
    ind2 = zeros(1,C(2));     
    
    for iT = 1:data.nT
        
        %write credit assignment before updating
        data_model.cat(iB,iT,1) = cat;
        
        %calculate net action values for each stage
        v1net = cat.*v1p1 + (1-cat).*v1p2;
        v2net = cat.*v2p1 + (1-cat).*v2p2;
        
        %calculate action probabilities
        q1 = it*v1net + pe*ind1; %Qvalue stage 1
        q2 = it*v2net + pe*ind2; %Qvalue stage 2
        ap1 = exp(q1 - logsumexp(q1,2)); %action probability stage 1
        ap2 = exp(q2 - logsumexp(q2,2)); %action probability stage 2
        
        %get rewards from data
        r(1) = data.r1(iB,iT);
        r(2) = data.r2(iB,iT);
        
        if strcmp(info.type,'behavioral fit') || strcmp(actiontype,'decision_given')
            c(1) = data.d1(iB,iT);
            c(2) = data.d2(iB,iT);
        elseif strcmp(info.type,'simulation')
            c(1) = fastrandsample(ap1);  % random choice
            c(2) = fastrandsample(ap2);  % random choice
            
            if strcmp(info.fbtype,'probs')
                r(1) = double(rand<data.r1(iB,iT,c(1)));
                r(2) = double(rand<data.r2(iB,iT,c(2)));
            elseif strcmp(info.fbtype,'value')
                r(1) = data.r1(iB,iT,c(1));
                r(2) = data.r2(iB,iT,c(2));
            end
            
        else
            warning('no type specified')
        end
        
        corr(1) = double(r(1)>0); %1 = correct, 0 = incorrect
        corr(2) = double(r(2)>0); %1 = correct, 0 = incorrect
        
        if any(isinf([log(ap1) log(ap2)]))
            gnu = 1;
        end
        
        %update log likelihood stage 1
        lik = lik + log(ap1(c(1))); %log likelihood
        %update log likelihood stage 2
        lik = lik + log(ap2(c(2))); %log likelihood
        
        
        %% model dynamics
        
        %calculate reward prediction error for policy 1
        RPEp1(1) = r(1)-v1p1(c(1)); RPEp1(2) = r(2)-v2p1(c(2));
        %calculate reward prediction error for policy 2
        RPEp2(1) = r(2)-v1p2(c(1)); RPEp2(2) = r(1)-v2p2(c(2));
        
        %calculate "hierarchical prediction error" for policy evaluation based on
        %surprise estimates
        surprise1 = sum(abs(RPEp1));
        surprise2 = sum(abs(RPEp2));
        HPE = surprise2 - surprise1;
        ca = ca + ar*HPE;
        cat = exp(ca)/(1+exp(ca));
        
        % update Qvalues
        % policy 1
        v1p1(c(1)) = v1p1(c(1)) + lr*RPEp1(1);
        v2p1(c(2)) = v2p1(c(2)) + lr*RPEp1(2);
        % policy 2
        v1p2(c(1)) = v1p2(c(1)) + lr*RPEp2(1);
        v2p2(c(2)) = v2p2(c(2)) + lr*RPEp2(2);
        
        % update perseveration
        ind1 = zeros(1,C(1)); ind2 = zeros(1,C(2));
        ind1(c(1)) = 1; ind2(c(2)) = 1;
                
        %forgetting
        v1p1(3-c(1)) = (1-lr)*v1p1(3-c(1));
        v2p1(3-c(2)) = (1-lr)*v2p1(3-c(2));
        v1p2(3-c(1)) = (1-lr)*v1p2(3-c(1));
        v2p2(3-c(2)) = (1-lr)*v2p2(3-c(2));
        
        %write Qvalues before updating
        data_model.v1p1_upd(iB,iT,:) = v1p1; data_model.v2p1_upd(iB,iT,:) = v2p1;
        data_model.v1p2_upd(iB,iT,:) = v1p2; data_model.v2p2_upd(iB,iT,:) = v2p2;
        data_model.cat_upd(iB,iT,1) = cat;
        
        %write simulated data structure (Qvalues before updating)
        data_model.c1(iB,iT,1) = c(1); data_model.c2(iB,iT,1) = c(2);
        data_model.r1(iB,iT,1) = r(1); data_model.r2(iB,iT,1) = r(2);
        data_model.corr1(iB,iT,1) = corr(1); data_model.corr2(iB,iT,1) = corr(2);
        data_model.RPEp1(iB,iT,:) = RPEp1; data_model.RPEp2(iB,iT,:) = RPEp2;
        data_model.HPE(iB,iT) = HPE;
        data_model.rpe_comp(iB,iT,:) = [[abs(RPEp2(1)) - abs(RPEp1(1))] [abs(RPEp2(2)) - abs(RPEp1(2))]];
        
        if any(isnan(v1p1)) || any(isnan(v2p1)) || any(isnan(v1p2)) || any(isnan(v2p2))
            warning('something in Qvalues is off')
        end
        
    end
end

if strcmp(info.type,'behavioral fit')
    output = lik;
elseif strcmp(info.type,'simulation')
    output  = data_model;
else
    error('no output specified')
end
% sumlik = mean(pred);


gnu = 1;
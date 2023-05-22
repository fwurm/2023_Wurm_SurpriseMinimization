function output = model_postCAT(x,data)

% Likelihood function for Q-learning on experiment postCAT v.3 with
% parameters (1)learning rate, (2)inverse temperature, (3)
% perseveration (4)start bias and (5)assignment rate
%
% USAGE: lik = model_postCAT(x,data)
%
% INPUTS:
%   x - parameters:
%       x(1) - learning rate
%       x(2) - inverse temperature
%       x(3) - perseveration
%       x(4) - start bias
%       x(5) - assignment rate
%   data - structure with the following fields
%          .d1 - [1 x nT] choices/decisions stage1
%          .d2 - [1 x nT] choices/decisions stage2
%          .r(1) - [1 x nT] rewards stage1
%          .r(2) - [1 x nT] rewards stage2
%          .N - number of trials (per block)
%
% OUTPUTS:
%   lik - log-likelihood
%
% Franz Wurm, December 2017
% adapted from Sam Gershman, June 2015

%% Pre settings stuff

info = [data.info];

if isfield(data,'d1') && isfield(data,'d2')
    actiontype = 'decision_given';
else
    actiontype = 'sample_decision';
end

if strcmp(info.type,'behavioral fit')
    C(1) = 2; % number of all options for stage 1
    C(2) = 2; % number of all options for stage 2
elseif strcmp(info.type,'simulation')
    C(1) = 2;
    C(2) = 2;
end


%get parameters
lr = x(1);
it = x(2);
pe = x(3);
ed = x(4); %evidence decay
dt = x(5); %decision threshold
%ca = x(6);
ca = 0.5;

%initialize
pred = zeros(1,2);
lik = 0;

%% Simulation
for iB = 1:data.nB
	
	e = 0; %initial evidence
	lastc = zeros(1,2); %choice on trial n-1
	lastr = zeros(1,2); %reward on trial n-1

	v1 = zeros(1,C(1));  % initial values stage 1
	v2 = zeros(1,C(2));  % initial values stage 1
	ind1 = zeros(1,C(1)); ind2 = zeros(1,C(2)); %indicator for perseveration



	for iT = 1:data.nT
    
		%write Qvalues before updating
		data_model.v1(iB,iT,:) = v1; data_model.v2(iB,iT,:) = v2;
    
		%calculate values and action values
		%stage 1
		q1 = it*v1 + pe*ind1; %Qvalue
		ap1 = exp(q1 - logsumexp(q1,2)); %action probability
		%stage 2
		q2 = it*v2 + pe*ind2; %Qvalue
		ap2 = exp(q2 - logsumexp(q2,2)); %action probability
    
		%get rewards from data
		r(1) = data.r1(iB,iT); 
		r(2) = data.r2(iB,iT);
    
		if strcmp(info.type,'behavioral fit') || strcmp(actiontype,'decision_given')
			c(1) = data.d1(iB,iT); 
			c(2) = data.d2(iB,iT);
		elseif strcmp(info.type,'simulation')
			c(1) = fastrandsample(ap1);  % random choice
			c(2) = fastrandsample(ap2);  % random choice
        
			if strcmp(info.fbtype,'prob')
				r(1) = double(rand<data.r1(iB,iT,c(1)));
				r(2) = double(rand<data.r2(iB,iT,c(2)));
			elseif strcmp(info.fbtype,'pay')
%             	r(1) = data.r1(iT,c(1));
%             	r(2) = data.r1(iT,c(2));
			elseif strcmp(info.fbtype,'relpay')
				r(1) = data.r1(iB,iT,c(1));
				r(2) = data.r2(iB,iT,c(2));
			end
        
        
%%% TODO: REWARD part
%         if strcmp(fbType,'prob')
%             r(1) = double(rand<R1(iT,c(1)));
%             r(2) = double(rand<R2(iT,c(2)));
%         elseif strcmp(fbType,'pay')
%             r(1) = R1(iT,c(1));
%             r(2) = R2(iT,c(2));
%         elseif strcmp(fbType,'relpay')
%             r(1) = 2*R1(iT,c(1)) - 1;
%             r(2) = 2*R2(iT,c(2)) - 1;
%         end
%%%               
		else
			warning('no type specified')
		end
    
		corr(1) = double(r(1)>0); %1 = correct, 0 = incorrect
		corr(2) = double(r(2)>0); %1 = correct, 0 = incorrect
    
		%update log likelihood stage 1    
		lik = lik + log(ap1(c(1))); %log likelihood    
		%update log likelihood stage 2   
		lik = lik + log(ap2(c(2))); %log likelihood
    
    
		%% model dynamics
    
		e = (1-ed)*e; %decay of evidence / forgetting
    
		if length(unique([c(1) c(2)] == lastc))==2 && length(unique([sign(r(1)) sign(r(2))] == lastr))==2
			if eq([c(1) c(2)],lastc)==eq([sign(r(1)) sign(r(2))],lastr)
				e = e + 1;
			else
				e = e - 1;
			end
		end
     
		%evidence-based reasoning
		if e > dt %H1
			ca = 1;
		elseif e < -dt %H2
			ca = 0;
		else %H0
			ca = 0.5;       
		end
    
    
    
		%calculate reward prediction error for stages
		rpe1(1) = r(1)-v1(c(1)); rpe1(2) = r(2)-v1(c(1));
		rpe2(1) = r(2)-v2(c(2)); rpe2(2) = r(1)-v2(c(2));
    
		% update values
		v1(c(1)) = v1(c(1)) + ca*lr*rpe1(1) + (1-ca)*lr*rpe1(2);
		v2(c(2)) = v2(c(2)) + ca*lr*rpe2(1) + (1-ca)*lr*rpe2(2);
        
        % forgetting
        v1(3-c(1)) = (1-lr)*v1(3-c(1));
        v2(3-c(2)) = (1-lr)*v2(3-c(2));
      
		% update perseveration
		ind1 = zeros(1,C(1)); ind2 = zeros(1,C(2));
		ind1(c(1)) = 1; ind2(c(2)) = 1;
    
		%Update last choice and reward
		lastc = c;
		lastr = [sign(r(1)) sign(r(2))];
		
		%write Qvalues after updating
		data_model.v1_upd(iB,iT,:) = v1; data_model.v2_upd(iB,iT,:) = v2;
    
		%write simulated data structure
		data_model.c1(iB,iT,1) = c(1); data_model.c2(iB,iT,1) = c(2);
		data_model.r1(iB,iT,1) = r(1); data_model.r2(iB,iT,1) = r(2);
		data_model.corr1(iB,iT,1) = corr(1); data_model.corr2(iB,iT,1) = corr(2);
		data_model.rpe1(iB,iT,:) = rpe1; data_model.rpe2(iB,iT,:) = rpe2;
		data_model.evidence(iB,iT,:) = e;
    
		if any(isnan(v1)) || any(isnan(v2))
			warning('something in Qvalues is off')
		end
    end
end


%% Output stuff
if strcmp(info.type,'behavioral fit')
    output = lik;
elseif strcmp(info.type,'simulation')
    output  = data_model;
else
    error('no output specified')
end
% sumlik = mean(pred);


gnu = 1;
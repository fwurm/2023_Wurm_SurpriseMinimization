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

info = [data.info];

C(1) = 2; % number of all options for stage 1
C(2) = 2; % number of all options for stage 2

%get parameters
it = x(1);

%initialize
lik = 0;


for iB = 1:data.nB
	v1 = zeros(1,C(1));  % initial values stage 1
	v2 = zeros(1,C(2));  % initial values stage 1

	for iT = 1:data.nT
    
		%write Qvalues before updating
		data_model.v1(iB,iT,:) = v1; data_model.v2(iB,iT,:) = v2;
    
		%calculate values and action values
		%stage 1
		q1 = it*v1; %Qvalue
		ap1 = exp(q1 - logsumexp(q1,2)); %action probability
		%stage 2
		q2 = it*v2; %Qvalue
		ap2 = exp(q2 - logsumexp(q2,2)); %action probability
    
		%get rewards from data
		r(1) = data.r1(iB,iT); 
		r(2) = data.r2(iB,iT);
    
		if strcmp(info.type,'behavioral fit')
			c(1) = data.d1(iB,iT); 
			c(2) = data.d2(iB,iT);
		elseif strcmp(info.type,'simulation')
			c(1) = fastrandsample(ap1);  % random choice
			c(2) = fastrandsample(ap2);  % random choice
        
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
    
		%calculate reward prediction error for stages
		rpe1(1) = r(1)-v1(c(1)); rpe1(2) = r(2)-v1(c(1));
		rpe2(1) = r(2)-v2(c(2)); rpe2(2) = r(1)-v2(c(2));
		
		%write Qvalues after updating
		data_model.v1_upd(iB,iT,:) = v1; data_model.v2_upd(iB,iT,:) = v2;
    
		%write simulated data structure (Qvalues before updating)
		data_model.c1(iB,iT,1) = c(1); data_model.c2(iB,iT,1) = c(2);
		data_model.r1(iB,iT,1) = r(1); data_model.r2(iB,iT,1) = r(2);
		data_model.corr1(iB,iT,1) = corr(1); data_model.corr2(iB,iT,1) = corr(2);
		data_model.rpe1(iB,iT,:) = rpe1; data_model.rpe2(iB,iT,:) = rpe2;
    
		if any(isnan(v1)) || any(isnan(v2))
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
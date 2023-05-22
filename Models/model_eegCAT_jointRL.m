function output = model_postCAT(x,data)

    % Likelihood function for Q-learning on experiment postCAT v.2 with
    % parameters (1)learning rate, (2)inverse temperature, (3)
    % perseveration and (4)credit assignment
    %
    % USAGE: lik = model_postCAT(x,data)
    %
    % INPUTS:
    %   x - parameters:
    %       x(1) - learning rate
    %       x(2) - inverse temperature
    %       x(3) - perseveration
    %   data - structure with the following fields
    %          .d1 - [1 x nT] choices/decisions stage1
    %          .d2 - [1 x nT] choices/decisions stage2
    %          .r1 - [1 x nT] rewards stage1
    %          .r2 - [1 x nT] rewards stage2
    %          .N - number of trials (per block)
    %
    % OUTPUTS:
    %   lik - log-likelihood
    %
    % Franz Wurm, December 2017
    % adapted from Sam Gershman, June 2015
    
    info = [data.info];
    
C = size(unique([data.d1]),1)*size(unique([data.d2]),1); % number of all options

lr = x(1);
b = x(2);
p = x(3);

lik = 0;

for iB = 1:data.nB

	v = zeros(1,C);  % initial values stage 1
	
	ind = zeros(1,C); %indicator for perseveration

	for iT = 1:data.nT
        
        %write Qvalues before updating
        data_model.v(iB,iT,:) = v; 
    
		c(1) = data.d1(iB,iT); c(2) = data.d2(iB,iT);
		r(1) = data.r1(iB,iT); r(2) = data.r2(iB,iT);
        rJ = r(1)+r(2);
    
		if c(1) == 1 && c(2) == 1
			cJ = 1;       
		elseif c(1) == 1 && c(2) == 2
			cJ = 2; 
		elseif c(1) == 2 && c(2) == 1
			cJ = 3; 
		elseif c(1) == 2 && c(2) == 2
			cJ = 4; 
		end  
    
		%update log likelihood
		q = b*v + p*ind; %Qvalue
		ap = exp(q - logsumexp(q,2)); %action probability
    
		lik_temp = log(ap(cJ));
		if isinf(log(ap))
			lik_temp = -(10^5);
		end
    
		lik = lik + lik_temp; %log likelihood
		pred(iT,1) = ap(cJ);

		%calculate reward prediction error
		rpe = rJ-v(cJ);
    
		% update values
		v(cJ) = v(cJ) + lr*rpe;
    
		% forgetting
		possib = [1:4];
		possib([cJ]) = [];
		v(possib) = (1-lr).*v(possib);
     
		% update perseveration
		ind = zeros(1,C);
		ind(cJ) = 1;
    
		if any(isinf(log(ap)))
			gnu = 1;;
        end
        
        %write Qvalues before updating
        data_model.v_upd(iB,iT,:) = v; 
        
        %write simulated data structure (Qvalues before updating)
        data_model.c(iB,iT,1) = c(1); data_model.c2(iB,iT,1) = c(2); data_model.cJ(iB,iT,1) = cJ;
        data_model.r1(iB,iT,1) = r(1); data_model.r2(iB,iT,1) = r(2); data_model.rJ(iB,iT,1) = rJ;
        data_model.corr1(iB,iT,1) = corr(1); data_model.corr2(iB,iT,1) = corr(2);
        data_model.rpeJ(iB,iT,:) = rpe; 
        
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
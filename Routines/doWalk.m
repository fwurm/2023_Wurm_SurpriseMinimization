function [payoff,dichot] = doWalk(options)

start = (options.upperBound-options.lowerBound).*rand(options.nB,options.nS) + options.lowerBound;

dichot = []; %draw from P(rew|correct)
reward_pre = start;

payoff = nan(options.nB,options.nS,options.nT);
payoff(:,:,1) = start;

for iB = 1:options.nB
    for iS = 1:options.nS
        for iT = 1:options.nT-1
            
            step = normrnd(options.walkMean,options.walkStd,[1 1]); %implies sampling ~normal(val_mean,val_Sd)
            
            reward_pre2(iB,iS) = reward_pre(iB,iS) + step; %add noise to walk
            
            while (reward_pre2(iB,iS) >= options.upperBound) || (reward_pre2(iB,iS) <= options.lowerBound) %check upper and lower bound
                step = normrnd(options.walkMean,options.walkStd,[1 1]);
                reward_pre2(iB,iS) = reward_pre(iB,iS) + step;
            end
            reward_pre(iB,iS) = reward_pre2(iB,iS); %new value get old value (random walk OHNE zurueklaufen)
            
            payoff(iB,iS,iT+1) = reward_pre2(iB,iS); %write new value
        end
    end
    
end

%dichotimization of reward
for iT = 1:options.nT
    draw = 100.*rand(options.nB,2);
    draw(draw<=payoff(:,:,iT)) = 1;
    draw(draw>=payoff(:,:,iT)) = 0;
    dichot(:,:,iT) = draw;
end
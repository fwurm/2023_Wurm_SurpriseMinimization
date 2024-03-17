function fw_modelcomp(dir,exp,models)

fprintf('\nModel comparison\n')


for i = 1:length(models)
    load(fullfile(dir.dir_fits,strcat('fit-',models{i},'-',exp,'.mat')),'results')
    fit(i) = results;
end

bms_results = mfit_bms(fit,0); %Laplace-informed priors
% bms_results = mfit_bms(fit,1); %BIC-informed priors


for i = 1:length(models)
    fprintf('   model %d: %s\n',i,models{i})
    fprintf('      alpha: %.2f\n',[bms_results.alpha(i)])
    fprintf('      expectation of p(m|y): %.2f\n',[bms_results.exp_r(i)])
    fprintf('      exceedence probability: %.2f\n',[bms_results.xp(i)])
    fprintf('      protected exceedence probability: %.2f\n',[bms_results.pxp(i)])
    fprintf('      BIC[mean(sd)]: %.2f(%.2f)\n',mean([fit(i).bic]),std([fit(i).bic]))
    fprintf('      BIC[sum]: %.2f\n',sum(sum([fit(i).bic])))
    fprintf('      AIC[mean(sd)]: %.2f(%.2f)\n',mean([fit(i).aic]),std([fit(i).aic]))
    fprintf('      AIC[sum]: %.2f\n',sum(sum([fit(i).aic])))
    fprintf('      -LL[sum]: %.2f\n',sum(sum([fit(i).loglik])))
    
    infomat(i,:) = [sum(sum([fit(i).loglik])) sum(sum([fit(i).bic])) sum(sum([fit(i).aic])) [bms_results.xp(i)]];
end

gnu = 1;

%% try model family analysis (results unclear)
% % options.families = {[1], [2:7] [8:17]};
% % options.families = {[1:6],[7:16]};
% % options.families = {[1 2],[3 4], [5 6]}; % yoked vs. free decay vs. counterfactual
% % options.families = {[1 3 4],[2 5 6]}; % normal vs. start
% % options.families = {[1 3 6 8],[2 4 7 9], [5 10]}; %correct vs incorrect vs random
% % options.families = {[1 2 6 7],[3 4 8 9], [5 10]}; %100 vs 90 vs random
% % [posterior, out] = VBA_groupBMC(bms_results.lme', options) ; %
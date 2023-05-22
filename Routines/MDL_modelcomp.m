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
end
function [OUT] = getBehavior_vp_eegCAT(IN,trialnums,set,pos)

%aggregate variables
valueP1 = [IN.value1_p1  IN.value2_p1];
valueP2 = [IN.value1_p2  IN.value2_p2];
difficultyP1 = [IN.diff1_p1  IN.diff2_p1];
difficultyP2 = [IN.diff1_p2  IN.diff2_p2];

%sort RPEs &Qupdated
nTrial = length(trialnums);
rpe_p1_sort = nan(nTrial,1);
rpe_p2_sort = nan(nTrial,1);
Q_p1_sort = nan(nTrial,1);
Q_p2_sort = nan(nTrial,1);
value_p1_sort = nan(nTrial,1);
value_p2_sort = nan(nTrial,1);
diff_p1_sort = nan(nTrial,1);
diff_p2_sort = nan(nTrial,1);
hpe_sort = nan(nTrial,1); %evidence
caw = nan(nTrial,1); %credit assignment weight
caw2 = nan(nTrial,1); %credit assignment weight

% caw_pre = reshape(repmat(IN.caw_,[1 2])',1,2*numel(IN.caw))';

for iT = 1:length(trialnums)
    rpe_p1_sort(iT,1) = IN.rpe_p1(trialnums(iT),set(iT));
    rpe_p2_sort(iT,1) = IN.rpe_p2(trialnums(iT),set(iT));
%     Q_p1_sort = [Q_p1_sort; QP1(trialnums(iT),seq(iT))];
%     Q_p2_sort = [Q_p2_sort; QP2(trialnums(iT),seq(iT))];
    value_p1_sort(iT,1) = valueP1(trialnums(iT),set(iT));
    value_p2_sort(iT,1) = valueP2(trialnums(iT),set(iT));
    diff_p1_sort(iT,1) = difficultyP1(trialnums(iT),set(iT));
    diff_p2_sort(iT,1) = difficultyP2(trialnums(iT),set(iT));
    hpe_sort(iT,1) = IN.hpe(trialnums(iT),set(iT));
    
    if pos(iT) == 2
        caw(iT,1) = IN.caw_new(trialnums(iT));
    elseif pos(iT) == 1
        cat = IN.caw_old(trialnums(iT));
        ca = log(cat/(1-cat)); %logit (inverse sigmoidal);
        surprise1 = abs(rpe_p1_sort(iT,1));
        surprise2 = abs(rpe_p2_sort(iT,1));
        evidence = surprise2 - surprise1;
        ca2 = ca + IN.param(4) * evidence;
        cat2 = exp(ca2)./(1+exp(ca2));
        caw(iT,1) = cat2;
    else
        error('something is off')
    end


%         caw(iT,1) = caw_pre(trialnums(iT));
end

% valence is the direction of the rpe
valence_p1 = sign(rpe_p1_sort);
valence_p2 = sign(rpe_p2_sort);

% surprise is the distance to zero
surprise_p1 = abs(rpe_p1_sort);
surprise_p2 = abs(rpe_p2_sort);

OUT.rpeP1 = rpe_p1_sort;
OUT.rpeP2 = rpe_p2_sort;
OUT.valenceP1 = valence_p1;
OUT.valenceP2 = valence_p2;
OUT.surpriseP1 = surprise_p1;
OUT.surpriseP2 = surprise_p2;
OUT.valueP1 = value_p1_sort;
OUT.valueP2 = value_p2_sort;
OUT.difficultyP1 = diff_p1_sort;
OUT.difficultyP2 = diff_p2_sort;
OUT.hpe = hpe_sort;
OUT.caw = caw;
OUT.set = set;
OUT.pos = pos;




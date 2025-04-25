% figure1_2_S1_plots_P1

% for participant P1

% plots the behavioral accuracy across sessions (Fig.1d)
% plots the percentage of trials felt, sorted by current amplitude (Fig.2a)
% plots the percentage of trials felt, sorted by current amplitude and offet timing (Fig.S1a)

% Isabelle Rosenthal 2025

% load data
load(['..\Data\stimContext_Timing_P1_preprocessedSpks.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots the behavioral accuracy across sessions (Fig.1d)

% color map
condsC = [93 187 71;
    22 95 58;
    0 0 0]/255;

% pull data from structure
ct = 0; % data block counter
tt = 0; ft = 0; % counter for trials, felt trials
ansMatch = []; % number of trials with correct answers
ansTotal = []; % total number of trials
runType = [];
for di = 1:size(dataStruct,1) % for each day
    for si = 1:size(dataStruct,2) % for each run
        if dataStruct{di,si}.runType ~=0 % not including baseline
            ct=ct+1;
            
            % get the trials where a sensation was elicited
            feltTrials = ~isnan(dataStruct{di,si}.trialIntensity);
            runType(ct) = dataStruct{di,si}.runType; % 1=realistic, 2=abstract
            
            senAns = dataStruct{di,si}.trialSenOrder;
            senOrder = nan(size(senAns));
            senOrder(strcmp(senAns,'Vision First')) = -1;
            senOrder(strcmp(senAns,'Stim First')) = 1;
            senOrder(strcmp(senAns,'Simultaneous')) = 0;
            realOrder = dataStruct{di,si}.trialVisRelativeToStim;
            realOrderBinary = realOrder;
            realOrderBinary(realOrder~=0) = realOrder(realOrder~=0)./abs(realOrder(realOrder~=0));
            ansMatch(ct) = sum(senOrder(feltTrials)==(realOrderBinary(feltTrials)));
            ansTotal(ct) = sum(feltTrials);
            
            for tr = 1:size(dataStruct{di,si}.isCatchTrials,2)
                if dataStruct{di,si}.isCatchTrials(tr)==0 % during non-catch trials
                    tt=tt+1;
                    isFelt(tt) = ~isnan(dataStruct{di,si}.trialIntensity(tr))+1; % 1=not felt, 2=felt
                    
                    if isFelt(tt)==2 % if the trial had an elicited sensation
                        ft = ft+1;
                        isAccurate(ft) = (senOrder(tr) == realOrderBinary(tr))+1; % 1=not right, 2=right
                        c_amps(ft) =  dataStruct{di,si}.trialAmps(tr);
                        c_time(ft) =  dataStruct{di,si}.trialVisRelativeToStim(tr);
                        c_cond(ft) =  dataStruct{di,si}.runType;
                    end
                end
            end
        end
    end
end

perAnsCorOverallDay = (ansMatch(runType==1) + ansMatch(runType==2))...
    ./(ansTotal(runType==1)+ansTotal(runType==2));
perAnsCorRealistic = ansMatch(runType==1)./ansTotal(runType==1);
perAnsCorAbstract = ansMatch(runType==2)./ansTotal(runType==2);

disp('P1 results: task accuracy')
disp(['realistic trials felt by session: mean= '...
    num2str(mean(ansTotal(runType==1))) ', std=' num2str(std(ansTotal(runType==1)))])
disp(['abstract trials felt by session: mean= '...
    num2str(mean(ansTotal(runType==2))) ', std=' num2str(std(ansTotal(runType==2)))])

% F-test: was there an effect of learning over time across sessions
disp('F-test on behavioral accuracy: was there an effect of learning over time across sessions?')
mdl = fitlm(1:9, perAnsCorOverallDay);
[p, F] = coefTest(mdl);
disp(['no effect: p=' num2str(p)])

% logistic regression test
c_X = [c_cond]';
[B_c, dev_c, stats_c] = mnrfit(c_X, isAccurate', 'model','nominal');
p_c = [stats_c.p(2)];
% the first term in B is the intercept
disp(['logistic regression test'])
disp(['significance of condition in giving correct answer: ' num2str(p_c)])

% plot the scatter plot of behavior accuracy, with line for the mean
figure;
hold on
plot([0.5 9.5],[1/3 1/3],'--','Color',[0.5 0.5 0.5],'LineWidth',2); % add chance level
fMh(1) =plot(perAnsCorOverallDay,'Color',condsC(3,:),'LineWidth',3);
fMh(2) = plot([1:size(dataStruct,1)]-0.1,perAnsCorRealistic,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor',condsC(1,:), 'MarkerSize',12);
fMh(3) =plot([1:size(dataStruct,1)]+0.1,perAnsCorAbstract,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor',condsC(2,:), 'MarkerSize',12);
ylim([0 1])
xlim([0.5 9.5])
xticks([1:9])
xlabel('session')
ylabel('% correct')
legend(fMh(1:3),{'mean'; 'realistic';'abstract'}, 'Location', 'northwest')
title('P1: performance by visual condition')

% plot accompanying histogram of behavior accuracy across days
figure;
hold on
fPh(1) = histogram(perAnsCorRealistic,'BinEdges',[0:0.05:1],'FaceColor',condsC(1,:),'FaceAlpha',0.7,'Normalization','probability');
fPh(2) = histogram(perAnsCorAbstract,'BinEdges',[0:0.05:1],'FaceColor',condsC(2,:),'FaceAlpha', 0.7,'Normalization','probability');
xlim([0 1])
ylim([0 0.6])
xlabel('% correct')
ylabel('% of data sessions')
legend(fPh,{'realistic';'abstract'}, 'Location', 'northwest')
title('P1: histogram performance by visual condition')

%% plots the percentage of trials felt, sorted by current amplitude (Fig.2a)
% AND plots the percentage of trials felt, sorted by current amplitude and offet timing (Fig.S1a)

rng(21); % set random seed for bootstrap
ct = [0 0 0]; % data block counter, by condition type
ampsUsed = [100 60 30 0]; % current amplitudes tested
timesUsed = [-300  -150     0   150   300 Inf]; % timing offsets tested for P1
allTimes = [-300:75:300]; % all times tested across both participants (for plotting)

% colormap
runHues = [0.549019608	0.549019608	0.549019608; % baseline
    0.08627451	0.37254902	0.22745098; % abstract
    0.364705882	0.733333333	0.278431373;];

isFeltAll = cell(2, 3,5);
trialsFelt = [];
trialsTotal = [];
trialsFeltBaseline = [];
trialsTotalBaseline = [];
tt = 0; ft = 0; bt= 0; % counter for trials, felt trials, baseline trials
for di = 1:size(dataStruct,1) % for each day
    for si = 1:size(dataStruct,2) % for each run
        rt = dataStruct{di,si}.runType;
        ct(rt+1)=ct(rt+1)+1;
        
        % get the trials where a sensation was elicited
        feltTrials = ~isnan(dataStruct{di,si}.trialIntensity); % trialsFelt = condition x sess x amp x times
        
        % count the number of trials felt by amplitude and timing
        if rt ~=0 % not including baseline
            for aI = 1:numel(ampsUsed)-1
                for tI = 1:numel(timesUsed)-1
                    tTrials = (dataStruct{di,si}.trialAmps==ampsUsed(aI)) & (dataStruct{di,si}.trialVisRelativeToStim==timesUsed(tI));
                    trialsFelt(rt,ct((rt+1)), aI,tI) = sum(tTrials & feltTrials);
                    trialsTotal(rt,ct((rt+1)), aI,tI) = sum(tTrials);
                    isFeltAll{rt,aI,tI} = [isFeltAll{rt,aI,tI}  feltTrials(tTrials)];
                end
            end
            
            % sort trials for statistical tests
            for tr = 1:size(dataStruct{di,si}.isCatchTrials,2)
                if dataStruct{di,si}.isCatchTrials(tr)==0 % during non-catch trials
                    tt=tt+1;
                    isFelt(tt) = ~isnan(dataStruct{di,si}.trialIntensity(tr))+1; % 1=not felt, 2=felt
                    f_amps(tt) =  dataStruct{di,si}.trialAmps(tr);
                    f_time(tt) =  dataStruct{di,si}.trialVisRelativeToStim(tr);
                    f_cond(tt) =  dataStruct{di,si}.runType;
                end
            end
            
        else % if baseline
            for aI = 1:numel(ampsUsed)-1
                tTrials = (dataStruct{di,si}.trialAmps==ampsUsed(aI));
                trialsFeltBaseline(ct(rt+1), aI) = sum(tTrials & feltTrials);
                trialsTotalBaseline(ct(rt+1), aI) = sum(tTrials);
            end
            
            % sort trials for statistical tests
            for tr = 1:size(dataStruct{di,si}.isCatchTrials,2)
                bt = bt+1;
                isFelt_B(bt) = ~isnan(dataStruct{di,si}.trialIntensity(tr))+1; % 1=not felt, 2=felt
                b_amps(bt) =  dataStruct{di,si}.trialAmps(tr);
            end
        end
    end
end
% get percentages
perTrialsFelt = squeeze(sum(trialsFelt,2)./sum(trialsTotal,2));
perTrialsFeltBaseline = sum(trialsFeltBaseline,1)./sum(trialsTotalBaseline,1);
perTrialsFelt(isnan(perTrialsFelt)) = -0.01;
perTrialsFeltBaseline(isnan(perTrialsFeltBaseline)) = -0.01;

perTrialsFeltByCondSummary = squeeze(sum(sum(trialsFelt,2),4)./sum(sum(trialsTotal,2),4)); % cond x amp

% bootstrap across the trials to generate 95% CIs
totalTrBaseline = sum(trialsTotalBaseline,1);
trFeltBaselineSum = sum(trialsFeltBaseline,1);
totalTr = squeeze(sum(trialsTotal,2));
trFeltSum = squeeze(sum(trialsFelt,2));
for aI = 1:numel(ampsUsed)-1
    trFeltVecBaseline = zeros(1,totalTrBaseline(aI));
    trFeltVecBaseline(1:trFeltBaselineSum(aI)) = 1;
    
    trFeltVecReal = zeros(1,sum(totalTr(1,aI,:)));
    trFeltVecReal(1:sum(trFeltSum(1,aI,:))) = 1;
    
    trFeltVecAbs = zeros(1,sum(totalTr(2,aI,:)));
    trFeltVecAbs(1:sum(trFeltSum(2,aI,:))) = 1;
    
    for ti = 1:5
        trFeltVecRealT{ti} = zeros(1,sum(totalTr(1,aI,ti)));
        trFeltVecRealT{ti}(1:sum(trFeltSum(1,aI,ti))) = 1;
        
        trFeltVecAbsT{ti} = zeros(1,sum(totalTr(2,aI,ti)));
        trFeltVecAbsT{ti}(1:sum(trFeltSum(2,aI,ti))) = 1;
    end
    for iter = 1:1000
        % boot baseline
        baseSamp = randsample(totalTrBaseline(1),totalTrBaseline(1),'true'); % same number of trials per amp
        realSamp = randsample(sum(totalTr(1,aI,:)),sum(totalTr(1,aI,:)),'true'); % same number of trials per amp
        absSamp = randsample(sum(totalTr(2,aI,:)),sum(totalTr(2,aI,:)),'true'); % same number of trials per amp
        
        % pull the trials
        perTrialsFeltBaselineBoot(aI, iter) = mean(trFeltVecBaseline(baseSamp));
        perTrialsFeltByCondSummaryBoot(1, aI, iter) = mean(trFeltVecReal(realSamp));
        perTrialsFeltByCondSummaryBoot(2, aI, iter) = mean(trFeltVecAbs(absSamp));
        
        % boot by individual times within real/abstract
        for ti = 1:5
            realSampT = randsample(totalTr(1,aI,ti),totalTr(1,aI,ti),'true'); % same number of trials per amp
            absSampT = randsample(totalTr(2,aI,ti),totalTr(2,aI,ti),'true'); % same number of trials per amp
            
            perTrialsFeltByCondSummaryBootT(1, aI,ti, iter) = mean(trFeltVecRealT{ti}(realSampT));
            perTrialsFeltByCondSummaryBootT(2, aI,ti, iter) = mean(trFeltVecAbsT{ti}(absSampT));
        end
        
    end
    perBaseBoot(aI,:) = [prctile(perTrialsFeltBaselineBoot(aI,:), 2.5) prctile(perTrialsFeltBaselineBoot(aI,:), 97.5)];
    summaryBoot(1,aI,:) = [prctile(perTrialsFeltByCondSummaryBoot(1,aI,:), 2.5) prctile(perTrialsFeltByCondSummaryBoot(1,aI,:), 97.5)];
    summaryBoot(2, aI,:) = [prctile(perTrialsFeltByCondSummaryBoot(2,aI,:), 2.5) prctile(perTrialsFeltByCondSummaryBoot(2,aI,:), 97.5)];
    
    for ti = 1:5
        summaryTBoot(1,aI,ti,:) = [prctile(perTrialsFeltByCondSummaryBootT(1,aI,ti,:), 2.5) prctile(perTrialsFeltByCondSummaryBootT(1,aI,ti,:), 97.5)];
        summaryTBoot(2, aI,ti,:) = [prctile(perTrialsFeltByCondSummaryBootT(2,aI,ti,:), 2.5) prctile(perTrialsFeltByCondSummaryBootT(2,aI,ti,:), 97.5)];
    end
end

disp('P1 results: rate of eliciting ICMS sensations')

% bar plot summary: the percentage of trials felt, sorted by current amplitude (Fig.2a)
figure('Position',[1000 150 300 450]);
for  aI = 1:numel(ampsUsed)-1
    barData = [perTrialsFeltBaseline(aI) perTrialsFeltByCondSummary(2,aI) perTrialsFeltByCondSummary(1,aI)];
    subplot(3,1,aI)
    b(aI) = bar([1 2 3],100*barData,...
        'FaceColor','flat');
    hold on
    for ind = 1:3
        b(aI).CData(ind,:) = runHues(ind,:);
    end
    lowerBar = 100*(barData - [perBaseBoot(aI,1) summaryBoot(2,aI,1) summaryBoot(1,aI,1)]);
    upperBar = 100*([perBaseBoot(aI,2) summaryBoot(2,aI,2) summaryBoot(1,aI,2)] - barData);
    errorbar([1 2 3],100*barData, lowerBar, upperBar,'.k');
    xticklabels({'baseline','abstract','realistic'})
    yticks([0:25:100])
    ylabel('percent trials felt')
    title([num2str(ampsUsed(aI)) 'A'])
    ylim([0 110])
    disp([num2str(ampsUsed(aI)) 'uA baseline mean % trials felt: ' num2str(100*perTrialsFeltBaseline(aI)), ...
        '% [' num2str(100*perBaseBoot(aI,1)) ', ' num2str(100*perBaseBoot(aI,2)) ']']);
    disp([num2str(ampsUsed(aI)) 'uA abstract mean % trials felt: ' num2str(100*perTrialsFeltByCondSummary(2,aI)), ...
        '% [' num2str(100*summaryBoot(2,aI,1)) ', ' num2str(100*summaryBoot(2,aI,2)) ']']);
    disp([num2str(ampsUsed(aI)) 'uA realistic mean % trials felt: ' num2str(100*perTrialsFeltByCondSummary(1,aI)), ...
        '% [' num2str(100*summaryBoot(1,aI,1)) ', ' num2str(100*summaryBoot(1,aI,2)) ']']);
end
sgtitle('P1')

% bar plot individual times: the percentage of trials felt, sorted by current amplitude and offet timing (Fig.S1a)
figure('Position',[1000 150 825 450]);
for  aI = 1:numel(ampsUsed)-1
    barData = [perTrialsFeltBaseline(aI) squeeze(perTrialsFelt(2,aI,:))' squeeze(perTrialsFelt(1,aI,:))'];
    subplot(3,1,aI)
    b(aI) = bar([0 1 2 3 4 5 6 7 8 9 10],100*[perTrialsFeltBaseline(aI) squeeze(perTrialsFelt(2,aI,:))' squeeze(perTrialsFelt(1,aI,:))'],...
        'FaceColor','flat');
    for ind = 1:11
        if ind == 1
            b(aI).CData(ind,:) = runHues(1,:);
        elseif ind<7
            b(aI).CData(ind,:) = runHues(2,:);
        else
            b(aI).CData(ind,:) = runHues(3,:);
        end
    end
    hold on
    errorbar([0 1 2 3 4 5 6 7 8 9 10],100*barData,...
        100*(barData - [perBaseBoot(aI,1) squeeze(summaryTBoot(2,aI,:,1))' squeeze(summaryTBoot(1,aI,:,1))']),...
        100*([perBaseBoot(aI,2) squeeze(summaryTBoot(2,aI,:,2))' squeeze(summaryTBoot(1,aI,:,2))'] - barData),'.k');
    xticks([0 1:0.5:10])
    xticklabels({'baseline',allTimes,allTimes})
    yticks([0:25:100])
    ylabel('percent trials felt')
    title([num2str(ampsUsed(aI)) 'A'])
    ylim([0 110])
end
sgtitle('P1')

% accompanying statistics for Fig.2a and Fig. S1a

% test baseline effect of amp on is felt
[BB, devb, statsb] = mnrfit(b_amps, isFelt_B', 'model','nominal');
p_b = [statsb.p(2)];
% the first term in B is the intercept
disp(['significance of amp in baseline whether a trial was felt: ' num2str(p_b)])

% logisitic regression test
% test realistic and abstract conditions, timing, amplitude
% effects on whether a sensation was felt

% combine X predictors together
f_X = [f_amps; f_time; f_cond]';
[B, dev, stats] = mnrfit(f_X, isFelt', 'model','nominal');
p_f = [stats.p(2:4)];
% bonferroni correct since we are combining/comparing categories
p_f = utils.bonf_holm(p_f);
% the first term in B is the intercept
disp(['significance of amp in whether a trial was felt: ' num2str(p_f(1))])
disp(['significance of timing in whether a trial was felt: ' num2str(p_f(2))])
disp(['significance of condition in whether a trial was felt: ' num2str(p_f(3))])

% logisitic regression test - within 60uA trials
% test baseline vs visual (combine realistic and abstract) conditions
% effects on whether a sensation was felt

% combine X predictors together
f_X = [ones(1,sum(b_amps==60))*1 ones(1,sum(f_amps==60))*2]';
feltV = [isFelt_B(b_amps==60)'; isFelt(f_amps==60)'];
[B, dev, stats] = mnrfit(f_X, feltV, 'model','nominal');
p_60 = [stats.p(2)];
% the first term in B is the intercept
disp(['in 60uA trials, significance of baseline vs visual condition in whether a trial was felt: ' num2str(p_60)])

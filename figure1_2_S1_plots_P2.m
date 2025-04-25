% figure1_2_S1_plots_P2

% for participant P2

% plots the behavioral accuracy across sessions (Fig.1d)
% plots the percentage of trials felt, sorted by current amplitude (Fig.2b)
% plots the percentage of trials felt, sorted by current amplitude and offet timing (Fig.S1b)

% Isabelle Rosenthal 2025

% load data
load(['..\Data\stimContext_Timing_P2_preprocessedSpks.mat'])

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
            feltTrials = ~strcmp(dataStruct{di,si}.trialSenOrder,'NaN');
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
                    isFelt(tt) = (~strcmp(dataStruct{di,si}.trialSenOrder(tr),'NaN'))+1; % 1=not felt, 2=felt
                    
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

disp('P2 results: task accuracy')
disp(['abstract accuracy over sessions: mean=' num2str(100*mean(perAnsCorAbstract)), ...
    '%, std=' num2str(100*std(perAnsCorAbstract)) '%']);
disp(['realistic accuracy over sessions: mean=' num2str(100*mean(perAnsCorRealistic)), ...
    '%, std=' num2str(100*std(perAnsCorRealistic)) '%']);

disp(['realistic trials felt by session: mean= '...
    num2str(mean(ansTotal(runType==1))) ', std=' num2str(std(ansTotal(runType==1)))])
disp(['abstract trials felt by session: mean= '...
    num2str(mean(ansTotal(runType==2))) ', std=' num2str(std(ansTotal(runType==2)))])

% F-test: was there an effect of learning over time across sessions
disp('F-test on behavioral accuracy: was there an effect of learning over time across sessions?')
mdl = fitlm(1:10, perAnsCorOverallDay);
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
plot([0.5 10.5],[2/5 2/5],'--','Color',[0.5 0.5 0.5],'LineWidth',2); % add chance level for stim first, vis first trials
plot([0.5 10.5],[1/5 1/5],'--','Color',[0.5 0.5 0.5],'LineWidth',2); % add chance level for simultaneous trials
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
title('P2: performance by visual condition')

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
title('P2: histogram performance by visual condition')

%% plots the percentage of trials felt, sorted by current amplitude (Fig.2b)
% AND plots the percentage of trials felt, sorted by current amplitude and offet timing (Fig.S1b)

rng(21); % set random seed for bootstrap
ct = [0 0 0]; % data block counter, by condition type
ampsUsed = [100 60 30 0]; % current amplitudes tested
timesUsed = [-300 -225 -150 -75    0 75  150 225  300 Inf]; % timing offsets tested for P2
allTimes = [-300:75:300]; % all times tested across both participants (for plotting)

% colormap
runHues = [
    0.08627451	0.37254902	0.22745098; % abstract
    0.364705882	0.733333333	0.278431373; % realistic
    ];

isFeltAll = cell(2, 3, size(timesUsed,2)-1);
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
        feltTrials =~strcmp(dataStruct{di,si}.trialSenOrder,'NaN'); % trialsFelt = condition x sess x amp x times
        
        % count the number of trials felt by amplitude and timing
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
                isFelt(tt) = (~strcmp(dataStruct{di,si}.trialSenOrder(tr),'NaN'))+1; % 1=not felt, 2=felt
                f_amps(tt) =  dataStruct{di,si}.trialAmps(tr);
                f_time(tt) =  dataStruct{di,si}.trialVisRelativeToStim(tr);
                f_cond(tt) =  dataStruct{di,si}.runType;
            end
        end
    end
end
% get percentages
perTrialsFelt = squeeze(sum(trialsFelt,2)./sum(trialsTotal,2));
perTrialsFelt(isnan(perTrialsFelt)) = -0.01;
perTrialsFeltByCondSummary = squeeze(sum(sum(trialsFelt,2),4)./sum(sum(trialsTotal,2),4)); % cond x amp

% bootstrap across the trials to generate 95% CIs
totalTr = squeeze(sum(trialsTotal,2));
trFeltSum = squeeze(sum(trialsFelt,2));
for aI = 1:numel(ampsUsed)-1
    
    trFeltVecReal = zeros(1,sum(totalTr(1,aI,:)));
    trFeltVecReal(1:sum(trFeltSum(1,aI,:))) = 1;
    
    trFeltVecAbs = zeros(1,sum(totalTr(2,aI,:)));
    trFeltVecAbs(1:sum(trFeltSum(2,aI,:))) = 1;
    
    for ti = 1:9
        trFeltVecRealT{ti} = zeros(1,sum(totalTr(1,aI,ti)));
        trFeltVecRealT{ti}(1:sum(trFeltSum(1,aI,ti))) = 1;
        
        trFeltVecAbsT{ti} = zeros(1,sum(totalTr(2,aI,ti)));
        trFeltVecAbsT{ti}(1:sum(trFeltSum(2,aI,ti))) = 1;
    end
    for iter = 1:1000
        % boot baseline
        realSamp = randsample(sum(totalTr(1,aI,:)),sum(totalTr(1,aI,:)),'true'); % same number of trials per amp
        absSamp = randsample(sum(totalTr(2,aI,:)),sum(totalTr(2,aI,:)),'true'); % same number of trials per amp
        
        % pull the trials
        perTrialsFeltByCondSummaryBoot(1, aI, iter) = mean(trFeltVecReal(realSamp));
        perTrialsFeltByCondSummaryBoot(2, aI, iter) = mean(trFeltVecAbs(absSamp));
        
        % boot by individual times within real/abstract
        for ti = 1:9
            realSampT = randsample(totalTr(1,aI,ti),totalTr(1,aI,ti),'true'); % same number of trials per amp
            absSampT = randsample(totalTr(2,aI,ti),totalTr(2,aI,ti),'true'); % same number of trials per amp
            
            perTrialsFeltByCondSummaryBootT(1, aI,ti, iter) = mean(trFeltVecRealT{ti}(realSampT));
            perTrialsFeltByCondSummaryBootT(2, aI,ti, iter) = mean(trFeltVecAbsT{ti}(absSampT));
        end
        
    end
    summaryBoot(1,aI,:) = [prctile(perTrialsFeltByCondSummaryBoot(1,aI,:), 2.5) prctile(perTrialsFeltByCondSummaryBoot(1,aI,:), 97.5)];
    summaryBoot(2, aI,:) = [prctile(perTrialsFeltByCondSummaryBoot(2,aI,:), 2.5) prctile(perTrialsFeltByCondSummaryBoot(2,aI,:), 97.5)];
    
    for ti = 1:9
        summaryTBoot(1,aI,ti,:) = [prctile(perTrialsFeltByCondSummaryBootT(1,aI,ti,:), 2.5) prctile(perTrialsFeltByCondSummaryBootT(1,aI,ti,:), 97.5)];
        summaryTBoot(2, aI,ti,:) = [prctile(perTrialsFeltByCondSummaryBootT(2,aI,ti,:), 2.5) prctile(perTrialsFeltByCondSummaryBootT(2,aI,ti,:), 97.5)];
    end
end

disp('P2 results: rate of eliciting ICMS sensations')

% bar plot summary: the percentage of trials felt, sorted by current amplitude (Fig.2a)
figure('Position',[1000 150 300 450]);
for  aI = 1:numel(ampsUsed)-1
    barData = [perTrialsFeltByCondSummary(2,aI) perTrialsFeltByCondSummary(1,aI)];
    subplot(3,1,aI)
    b(aI) = bar([0 1 2],[0 100*barData],...
        'FaceColor','flat');
    hold on
    for ind = 2:3
        b(aI).CData(ind,:) = runHues(ind-1,:);
    end
    lowerBar = 100*(barData - [summaryBoot(2,aI,1) summaryBoot(1,aI,1)]);
    upperBar = 100*([summaryBoot(2,aI,2) summaryBoot(1,aI,2)] - barData);
    errorbar([1 2],100*barData, lowerBar, upperBar,'.k');
    xticks([1:2])
    yticks([0:25:100])
    xticklabels({'abstract','realistic'})
    ylabel('percent trials felt')
    title([num2str(ampsUsed(aI)) 'A'])
    ylim([0 110])
    disp([num2str(ampsUsed(aI)) 'uA abstract mean % trials felt: ' num2str(100*perTrialsFeltByCondSummary(2,aI)), ...
        '% [' num2str(100*summaryBoot(2,aI,1)) ', ' num2str(100*summaryBoot(2,aI,2)) ']']);
    disp([num2str(ampsUsed(aI)) 'uA realistic mean % trials felt: ' num2str(100*perTrialsFeltByCondSummary(1,aI)), ...
        '% [' num2str(100*summaryBoot(1,aI,1)) ', ' num2str(100*summaryBoot(1,aI,2)) ']']);
end
sgtitle('P2')

figure('Position',[1000 150 825 450]);
for  aI = 1:numel(ampsUsed)-1
    barData = [squeeze(perTrialsFelt(2,aI,:))' squeeze(perTrialsFelt(1,aI,:))'];
    subplot(3,1,aI)
    b(aI) = bar([0:18],100*[0 squeeze(perTrialsFelt(2,aI,:))' squeeze(perTrialsFelt(1,aI,:))'],...
        'FaceColor','flat');
    for ind = 2:19
        if ind<11
            b(aI).CData(ind,:) = runHues(1,:);
        else
            b(aI).CData(ind,:) = runHues(2,:);
        end
    end
    hold on
    errorbar([1:18],100*barData,...
        100*(barData - [squeeze(summaryTBoot(2,aI,:,1))' squeeze(summaryTBoot(1,aI,:,1))']),...
        100*([squeeze(summaryTBoot(2,aI,:,2))' squeeze(summaryTBoot(1,aI,:,2))'] - barData),'.k');
    xticks([1:18])
    yticks([0:25:100])
    xticklabels({[timesUsed(1:end-1),timesUsed(1:end-1)]})
    ylabel('percent trials felt')
    title([num2str(ampsUsed(aI)) 'A'])
    ylim([0 110])
end
sgtitle('P2')

% accompanying statistics for Fig.2a and Fig. S1a

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
% figureS2ab_plots_P2

% for participant P2

% plots the features of eye tracking data (Fig.S2a)
% plots the number of tuned eye features (Fig.S2b)

% load data
load(['..\Data\stimContext_Timing_P2_preprocessedSpks.mat'])

%% format and select valid eye tracking trials to analyse

useEyeThresh = 0.6; % percent of time that eye tracker must be valid to use a trial
cA = 0; cR = 0; % abstract and realistic counters
absRuns = []; realRuns = []; % index of visual conditions
perA = []; perR = []; % percentage of usable time
clear catchDataUseA catchDataUseR % eye tracking data to use
for sesh = 1:size(eyeData,1) % for every set collected
    
    % get the session types
    realRuns(sesh) = find(cellfun(@(x) x.runType, dataStruct(sesh,:))==1);
    absRuns(sesh) = find(cellfun(@(x) x.runType, dataStruct(sesh,:))==2);
    % select the catch trials
    realCatch{sesh,:} = find(dataStruct{sesh,realRuns(sesh)}.isCatchTrials==1);
    absCatch{sesh,:} = find(dataStruct{sesh,absRuns(sesh)}.isCatchTrials==1);
    
    % process realistic condition
    catchDataR = squeeze(eyeData(sesh,realRuns(sesh),realCatch{sesh,:},:));
    startVisR = dataStruct{sesh,realRuns(sesh)}.visStart(realCatch{sesh,:});
    stimPhaseStartR = dataStruct{sesh,realRuns(sesh)}.stimPhaseStartTime(realCatch{sesh,:});
    % pull catch trials
    for ci = 1:size(catchDataR,1) %for every catch trial
        % get values from eye tracker
        eyeIsValid = double(catchDataR(ci).eyeIsValid_bn); % eye tracking valid
        leftEyeOpen = double(catchDataR(ci).leftEyeOpen_bn);  % left eye is open
        rightEyeOpen = double(catchDataR(ci).rightEyeOpen_bn); % right eye is open
        
        % pull 0.5s of ITI before vis start, and then all 3 segments of
        % animation
        phBins = [];
        [~, phBins(1)] = min(abs((startVisR(ci)-0.5) -  catchDataR(ci).binEdges(:,1)));
        [~, phBins(2)] = min(abs((startVisR(ci)) -  catchDataR(ci).binEdges(:,1)));
        [~, phBins(3)] = min(abs((startVisR(ci)+0.5) -  catchDataR(ci).binEdges(:,1)));
        [~, phBins(4)] = min(abs((startVisR(ci)+1) -  catchDataR(ci).binEdges(:,1)));
        [~, phBins(5)] = min(abs((startVisR(ci)+1.5) -  catchDataR(ci).binEdges(:,1)));
        
        % calculate the amount of time that is usable
        for ph = 1:4
            seg = eyeIsValid(phBins(ph):phBins(ph+1)) & (leftEyeOpen(phBins(ph):phBins(ph+1)) | rightEyeOpen(phBins(ph):phBins(ph+1)));
            perR(sesh,ci,ph) = mean(seg);
        end
        
        % if usable, put in a new structure
        if all(perR(sesh,ci,:)>= useEyeThresh)
            cR = cR+1;
            catchDataR(ci).startVis = startVisR(ci);
            catchDataUseR(cR) = catchDataR(ci);
            catchDataUseR_sess(cR) = sesh; % set the data came from
            catchDataUseR_stimPhaseStart(cR) = stimPhaseStartR(ci);  % start time
        end
    end
    
    % process abstract condition
    catchDataA = squeeze(eyeData(sesh,absRuns(sesh),absCatch{sesh,:},:));
    startVisA = dataStruct{sesh,absRuns(sesh)}.visStart(absCatch{sesh,:});
    stimPhaseStartA = dataStruct{sesh,absRuns(sesh)}.stimPhaseStartTime(absCatch{sesh,:});
    % pull catch trials
    for ci = 1:size(catchDataA,1) %for every catch trial
        % get values from eye tracker
        eyeIsValid = double(catchDataA(ci).eyeIsValid_bn); % eye tracking valid
        leftEyeOpen = double(catchDataA(ci).leftEyeOpen_bn); % left eye is open
        rightEyeOpen = double(catchDataA(ci).rightEyeOpen_bn); % right eye is open
        
        % pull 0.5s of ITI before vis start, and then all 3 segments of
        % animation
        phBins = [];
        [~, phBins(1)] = min(abs((startVisA(ci)-0.5) -  catchDataA(ci).binEdges(:,1)));
        [~, phBins(2)] = min(abs((startVisA(ci)) -  catchDataA(ci).binEdges(:,1)));
        [~, phBins(3)] = min(abs((startVisA(ci)+0.5) -  catchDataA(ci).binEdges(:,1)));
        [~, phBins(4)] = min(abs((startVisA(ci)+1) -  catchDataA(ci).binEdges(:,1)));
        [~, phBins(5)] = min(abs((startVisA(ci)+1.5) -  catchDataA(ci).binEdges(:,1)));
        
        % calculate the amount of time that is usable
        for ph = 1:4
            seg = eyeIsValid(phBins(ph):phBins(ph+1)) & (leftEyeOpen(phBins(ph):phBins(ph+1)) | rightEyeOpen(phBins(ph):phBins(ph+1)));
            perA(sesh,ci,ph) = mean(seg);
        end
        
        % if usable, put in a new structure
        if all(perA(sesh,ci,:)>= useEyeThresh)
            cA = cA+1;
            catchDataA(ci).startVis = startVisA(ci);
            catchDataUseA(cA) = catchDataA(ci);
            catchDataUseA_sess(cA) = sesh; % set the data came from
            catchDataUseA_stimPhaseStart(cA) = stimPhaseStartA(ci); % start time
        end
    end
end

disp(['threshold: ' num2str(useEyeThresh)]);
disp(['P2 abstract trials usable: ' num2str(cA) '/36 = ' num2str(100*cA/36) '%'])
disp(['P2 realistic trials usable: ' num2str(cR) '/36 = ' num2str(100*cR/36) '%'])

% normalize all channels within each run by dividing by the baseline ITI mean across trials
% normalize abstract trials
uSess = unique(catchDataUseA_sess);
for sI = 1:numel(uSess) % for every session
    sess = uSess(sI);
    trThisSess = find(catchDataUseA_sess==sess);
    chBaselines = [];
    for trial = 1:size(trThisSess,2) % get the ITI bins for each trial
        tmp = (catchDataUseA_stimPhaseStart(trThisSess(trial))...
            - catchDataUseA(trThisSess(trial)).binEdges(:,1)); % end at the last bin to be completely in the ITI
        tmp = tmp(tmp>=0); % only consider bins which began before the start of the stim phase
        [~, endITIb] = min(tmp);
        ITIb = [1 endITIb]; % start at the beginning of the trial
        chBaselines(trial,1:3) = nanmean(catchDataUseA(trThisSess(trial)).eyeOrigin_bn(ITIb(1):ITIb(2),:)); % average over baseline period
        chBaselines(trial,4:6) = nanmean(catchDataUseA(trThisSess(trial)).eyeVec_bn(ITIb(1):ITIb(2),:)); % average over baseline period
    end
    % average across all trials in a run
    chBaseFR_A(sI,:) = nanmean(chBaselines,1);
    chBaseFR_A(sI,(chBaseFR_A(sI,:) == 0)) = 0.000001; % avoid divide by zero errors
    for trial = 1:size(trThisSess,2)
        catchDataUseA(trThisSess(trial)).eyeOrigin_bnN = catchDataUseA(trThisSess(trial)).eyeOrigin_bn./chBaseFR_A(sI,1:3);
        catchDataUseA(trThisSess(trial)).eyeVec_bnN = catchDataUseA(trThisSess(trial)).eyeVec_bn./chBaseFR_A(sI,4:6);
    end
end

% normalize realistic trials
uSess = unique(catchDataUseR_sess);
for sI = 1:numel(uSess)  % for every session
    sess = uSess(sI);
    trThisSess = find(catchDataUseR_sess==sess);
    chBaselines = [];
    for trial = 1:size(trThisSess,2) % get the ITI bins for each trial
        tmp = (catchDataUseR_stimPhaseStart(trThisSess(trial))...
            - catchDataUseR(trThisSess(trial)).binEdges(:,1));% end at the last bin to be completely in the ITI
        tmp = tmp(tmp>=0); % only consider bins which began before the start of the stim phase
        [~, endITIb] = min(tmp);
        ITIb = [1 endITIb]; % start at the beginning of the trial
        chBaselines(trial,1:3) = nanmean(catchDataUseR(trThisSess(trial)).eyeOrigin_bn(ITIb(1):ITIb(2),:)); % average over baseline period
        chBaselines(trial,4:6) = nanmean(catchDataUseR(trThisSess(trial)).eyeVec_bn(ITIb(1):ITIb(2),:)); % average over baseline period
    end
    % average across all trials in a run
    chBaseFR_R(sI,:) = nanmean(chBaselines,1);
    chBaseFR_R(sI,(chBaseFR_R(sI,:) == 0)) = 0.000001; % avoid divide by zero errors
    for trial = 1:size(trThisSess,2)
        catchDataUseR(trThisSess(trial)).eyeOrigin_bnN = catchDataUseR(trThisSess(trial)).eyeOrigin_bn./chBaseFR_R(sI,1:3);
        catchDataUseR(trThisSess(trial)).eyeVec_bnN = catchDataUseR(trThisSess(trial)).eyeVec_bn./chBaseFR_R(sI,4:6);
    end
end
%% plots the features of eye tracking data (Fig.S2a)

taskTimings = [0 0.5 1 1.5]; % in s, [start visual (down), start touch, end touch, end visual (up)]

doSmooth = 1; % to perform smoothing on firing rates
smoothFactor = 5; % smoothing factor
alpha = 0.4; % transparency of plotted lines

% colormap
condHues = [0.364705882	0.733333333	0.278431373 % realistic
    0.08627451	0.37254902	0.22745098]; % abstract
sideColor = [0.7 0.7 0.7];
stimColor = [0.5 0.5 0.5];

binTs = [-0.7 3]; % time bins (s) relative to when animation starts

% get data by condition
eyeSeg = {};
for cati = 1:2 % for each visual condition (1=realistic, 2=abstract)
    if cati == 1
        cdata = catchDataUseR;
    else
        cdata = catchDataUseA;
    end
    % extract the bins for each trial
    binWin = [];
    for tr = 1: size(cdata,2)
        winTs = cdata(tr).startVis+[binTs(1) binTs(end)];
        
        [m1, b1] = min(abs(cdata(tr).binEdges(:,1)-winTs(1)));
        [m2, b2] = min(abs(cdata(tr).binEdges(:,1)-winTs(2)));
        assert(all([m1 m2]<0.0001)); % bins should be almost identical
        binWin(tr,:) = [b1 b2];
        assert(std(binWin(:,1)-binWin(:,2))==0); % check all bin windows have the same length
        eyeSeg{cati}(tr,:,:) = [cdata(tr).eyeOrigin_bnN(binWin(tr,1):binWin(tr,2),:)'
            cdata(tr).eyeVec_bnN(binWin(tr,1):binWin(tr,2),:)']'; % tr x timebin x ch
    end
end

for ch = 1:6 % for each channel to plot
    
    figure('Position', [200 100 480 425], 'Color',[1 1 1]);
    ylim([-1 2.5]);yy = ylim();
    yp = [yy(1) yy(1) yy(2) yy(2)];
    patch([taskTimings(2:3) taskTimings(3:-1:2)]-0.5,yp,stimColor, 'FaceAlpha',0.5, 'LineStyle','none')
    hold on
    patch([taskTimings(1:2) taskTimings(2:-1:1)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    patch([taskTimings(3:4) taskTimings(4:-1:3)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    
    % get data and plot
    for cati = 2:-1:1
        
        eyeCh = squeeze(eyeSeg{cati}(:,:,ch));
        
        % smooth data
        if doSmooth
            eyeCh = sgolayfilt(eyeCh,1,smoothFactor,[],2);
        end
        
        % extract data, removing nans
        a = 1:size(eyeCh,2);
        b = nanmean(eyeCh,1);
        c = 1.96*(nanstd(eyeCh,[],1)/sqrt(size(eyeCh,1)));
        keepIndex = ~isnan(b) & ~isnan(c);
        a = a(keepIndex);
        b = b(keepIndex);
        c = c(keepIndex);
        
        % plot
        utils.plotsem_alpha_ul([binTs(1):0.05:binTs(end)]-0.5,b,c,c,condHues(cati,:),condHues(cati,:),alpha, '-',5);
    end
    xlim([-1.2 2.5])
    
    title(['P2: eye feature ' num2str(ch)]);
end
%% run the tuning analysis and bootstrap

iter = 1000; % number of bootstrap iterations to run
runNames = {'realistic','abstract'};

% window to calculate tuning
binWidth = 0.1;
binTs = [-0.7:binWidth:3]; % (s) relative to when vis cue starts in unity
numBins = size(binTs,2)-1;
% window for rest calculation
restTs = [-1.75 -0.75]; % (s) relative to when vis cue starts in unity

% get data by condition
for cati = 1:2 % for each visual condition (1=realistic, 2=abstract)
    if cati == 1
        cdata = catchDataUseR;
    else
        cdata = catchDataUseA;
    end
    
    % extract the bins for each trial
    binWin = [];
    chBaselines = [];  % baseline firing rates by trial
    ITIWin = []; eyeSeg = [];
    for tr = 1: size(cdata,2)
        baselineBinTs = cdata(tr).startVis+restTs;
        
        % pull baseline bins
        [m1, b1] = min(abs(cdata(tr).binEdges(:,1)-baselineBinTs(1)));
        [m2, b2] = min(abs(cdata(tr).binEdges(:,1)-baselineBinTs(2)));
        assert(all([m1 m2]<0.0001)); % bins should be almost identical
        ITIWin(tr,:) = [b1 b2];
        chBaselines(:,tr) = nanmean([cdata(tr).eyeOrigin_bnN(ITIWin(tr,1):ITIWin(tr,2),:)'
            cdata(tr).eyeVec_bnN(ITIWin(tr,1):ITIWin(tr,2),:)']'); % ch x tr
        
        % pull tuning window bins
        for bn =  1:numBins
            winTs = cdata(tr).startVis+[binTs(bn) binTs(bn+1)];
            
            [m1, b1] = min(abs(cdata(tr).binEdges(:,1)-winTs(1)));
            [m2, b2] = min(abs(cdata(tr).binEdges(:,1)-winTs(2)));
            assert(all([m1 m2]<0.0001)); % bins should be almost identical
            binWin(tr,:) = [b1 b2];
            % average across bins within each window
            eyeSeg(:,tr,bn) = nanmean([cdata(tr).eyeOrigin_bnN(binWin(tr,1):binWin(tr,2),:)'
                cdata(tr).eyeVec_bnN(binWin(tr,1):binWin(tr,2),:)']'); % ch x trial x class (timebin)
        end
    end
    
    % get baselines averaged across trials
    spk_baseline = nanmean(chBaselines,2);
    
    % one hot code the labels
    mod_1HC = zeros(size(eyeSeg,2)*(numBins), numBins);
    ct = 0; segTall = [];
    for tr = 1: size(eyeSeg,2) % for each trial
        for li = 1: (size(eyeSeg,3)) %for each label not baseline
            ct = ct+1;
            mod_1HC(ct, li)=1;
            segTall(:,ct) = eyeSeg(:,tr,li);
        end
    end
    
    % append the ITI FR, duplicated to match the number of repetitions
    modBase_1HC = [mod_1HC; zeros(size(eyeSeg,2),numBins)];
    for ch = 1:6 %for each channel
        % append baseline FR
        segTallBase = [squeeze(segTall(ch,:))'; ones(size(eyeSeg,2), 1)...
            *spk_baseline(ch)];
        mdl = fitlm(modBase_1HC,segTallBase);
        pLM_mod_byBin(ch,cati,:) = mdl.Coefficients.pValue(2:end); % not getting offset beta
    end
    % correct for multiple comparisons
    for te = 1:size(pLM_mod_byBin, 3) % for each variable
        pLM_mod_byBin(:, cati,te) = bonf_holm(pLM_mod_byBin(:, cati,te));
    end
    
    % bootstrap
    disp(['beginning ' runNames{cati} ' tuning bootstrap'])
    for bt = 1:iter
        if rem(bt,200)==0
            disp(['on iteration ' num2str(bt) ' of ' num2str(iter)]);
        end
        % samples with replacement
        trToUseBoot = datasample(1:size(eyeSeg,2),size(eyeSeg,2)); % samples with replacement
        
        % one hot code the labels
        mod_1HC = zeros(size(eyeSeg,2)*(numBins), numBins);
        ct = 0; segTall = [];
        for tr = 1: size(eyeSeg,2) % for each trial
            for li = 1: (size(eyeSeg,3)) %for each label not baseline
                ct = ct+1;
                mod_1HC(ct, li)=1;
                segTall(:,ct) = eyeSeg(:,trToUseBoot(tr),li);
            end
        end
        
        % append the ITI FR, duplicated to match the number of repetitions
        modBase_1HC = [mod_1HC; zeros(size(eyeSeg,2),numBins)];
        for ch = 1:6 %for each channel
            % recompute baseline
            spk_baselineB = nanmean(chBaselines(ch,trToUseBoot),2);
            
            % append baseline FR
            segTallBase = [squeeze(segTall(ch,:))'; ones(size(eyeSeg,2), 1) *spk_baselineB];
            mdl = fitlm(modBase_1HC,segTallBase);
            pLM_mod_byBin_boot(ch,cati,:,bt) = mdl.Coefficients.pValue(2:end); % not getting offset beta
        end
        % correct for multiple comparisons
        for te = 1:size(pLM_mod_byBin, 3) % for each variable
            pLM_mod_byBin_boot(:, cati,te,bt) = bonf_holm(pLM_mod_byBin_boot(:, cati,te,bt));
        end
        
    end
end

% compute CIs for the boot
% get bootstrap values
isSigBoot = squeeze(sum(pLM_mod_byBin_boot<0.05,1));
isTunedUpperCI = prctile(isSigBoot,97.5,3);
isTunedLowerCI = prctile(isSigBoot,2.5,3);

%% plots the number of tuned eye features (Fig.S2b)

alpha = 0.4; % transparency of plot
taskTimings = [0 0.5 1 1.5]; % in s, [start visual (down), start touch, end touch, end visual (up)]

% colormap
condHues = [0.364705882	0.733333333	0.278431373 % realistic
    0.08627451	0.37254902	0.22745098]; % abstract

sideColor = [0.7 0.7 0.7];
stimColor = [0.5 0.5 0.5];

% get the times associated with each bin
binWinT = [];
for ind = 1:numBins
    binWinT(ind,:) = [binTs(ind) binTs(ind+1)];
    binCenters(ind) = mean([binTs(ind) binTs(ind+1)]);
end
binCenterTouch = binCenters-0.5; % to align to movement of touch onset, not vis start

% plot tuning for each condition
figure('Position',[200 100 950 375])
for cati = 1:2
    subplot(1,2,cati)
    
    % extract values to plot where p-values are significant
    thisP = squeeze(pLM_mod_byBin(:, cati,:));
    sumP = sum(thisP<0.05);
    lowerCI = squeeze(isTunedLowerCI(cati,:));
    upperCI = squeeze(isTunedUpperCI(cati,:));
    upperP = upperCI-sumP;
    lowerP = sumP-lowerCI;
    
    % patch for when visual animations took place
    ylim([0 7]);yy = ylim();
    yp = [yy(1) yy(1) yy(2) yy(2)];
    patch([taskTimings(2:3) taskTimings(3:-1:2)]-0.5,yp,stimColor, 'FaceAlpha',0.5, 'LineStyle','none')
    hold on
    patch([taskTimings(1:2) taskTimings(2:-1:1)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    patch([taskTimings(3:4) taskTimings(4:-1:3)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    
    [h0, h1(cati)] = utils.plotsem_alpha_ul(binCenterTouch,sumP,...
        lowerP, upperP,condHues(cati,:),condHues(cati,:),alpha, '-', 3);
    set(h0,'edgec',condHues(cati,:));
    plot(binCenterTouch,sumP,'o', 'LineWidth',3,...
        'Color',condHues(cati,:), 'MarkerFaceColor',condHues(cati,:),'MarkerSize',4);
    
    title([runNames{cati}]);
    ylim([0 7])
    xlim([binCenterTouch(1) binCenterTouch(end)])
    xticks([-1:0.5:2])
    
end
sgtitle(['P2: EYE linear regression tuned channels relative to baseline'])




% figure5abc_plots_P2

% for participant P2

% plots the number of tuned channels (Fig.5a)
% plots example channel activity (Fig.5b)
% plots the venn diagram of overlap of tuned channels between visual
% conditions (Fig.5c)

% Isabelle Rosenthal 2025

% if you want to load in the bootstrap analysis with 1000 iterations
% instead of rerunning it yourself
toLoad = 1; % 1=load in prerun analysis, 0=run it yourself
% if you run yourself, bootstrapped values may be slightly different due to
% change in the random seed

% load data
load(['..\Data\stimContext_Timing_P2_preprocessedSpks.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the tuning analysis and bootstrap
% OR load in the prerun file

iter = 1000; % number of bootstrap iterations to run
runNames = {'realistic','abstract'};

% window to calculate tuning
binWidth = 0.1;
binTs = [-0.7:binWidth:3]; % (s) relative to when vis cue starts in unity
numBins = size(binTs,2)-1;
% window for rest calculation
restTs = [-1.75 -0.75]; % (s) relative to when vis cue starts in unity

% extract normalized firing rates
dataN =  cellfun(@(x) x.dataN, dataStruct, 'UniformOutput',0);

% pull the correct blocks for each condition
runCats{1} = find(cellfun(@(x) x.runType, dataStruct)==1);
runCats{2} = find(cellfun(@(x) x.runType, dataStruct)==2);

if toLoad == 1
    % if you're loading in former results
    load(['..\Data\stimContextTiming_P2_tunedCatch_0.1bn_1000BOOTResults.mat']);
    disp('results loaded successfully')
    disp('**********************************')
else
    disp('Note that bootstrapped values may vary slightly due to randomness in the distribution.')
    disp('Set a random seed to ensure replicability.')
    
    % format into trial data and labels
    pLM_mod_byBin = []; % p values for tuning analysis
    pLM_mod_byBin_boot = []; % bootstrapped p values
    segments = []; % data by trial
    spk_baseline = []; % baseline firing rates averaged across trials
    chBaselines = []; % baseline firing rates by trial
    for cati = 1:2 % for each visual condition (1=realistic, 2=abstract)
        data = (dataN(runCats{cati})); data = cat(1,data{:});
        dataAmps = cellfun(@(x) x.trialAmps, dataStruct(runCats{cati}), 'UniformOutput',0);
        dataAmps = cat(2,dataAmps{:});
        dataBinsT = cellfun(@(x) x.binsT, dataStruct(runCats{cati}), 'UniformOutput',0);
        dataBinsT = cat(2,dataBinsT{:});
        uVisStartT = cellfun(@(x) x.visStart, dataStruct(runCats{cati}), 'UniformOutput',0);
        uVisStartT = cat(2,uVisStartT{:});
        
        % extract catch trials only
        tmpInd = (dataAmps==0);
        data = data(tmpInd,:,:);
        dataBinsT = dataBinsT(tmpInd);
        uVisStartT = uVisStartT(tmpInd);
        
        % combine both S1 arrays
        dataS1 = [data(:,:,1) data(:,:,2)];
        
        % extract the bins for each trial
        for tr = 1: size(dataS1,1)
            
            % get the times to use for baseline and main analysis
            % trials are aligned to visStart in preprocessing        
            baselineBinTs = uVisStartT(tr)+restTs;
            tuneBinTs = uVisStartT(tr)+binTs;
            
            spk =  squeeze(dataS1(tr,:)); % data for this trial
             
            % pull baseline bins
            [m1, b1] = min(abs(dataBinsT{tr}-baselineBinTs(1)));
            [m2, b2] = min(abs(dataBinsT{tr}-baselineBinTs(2)));
            assert(all([m1 m2]<0.0001)); % bins should be almost identical
            ITIb = [b1 b2];
            chBaselines(tr,:, cati) = cell2mat(cellfun(@(x) nanmean(x(ITIb(1):ITIb(2))), spk, 'UniformOutput', false));
            
            % pull tuning window bins
            for ind = 1:numBins
                binWinT(ind,:) = [tuneBinTs(ind) tuneBinTs(ind+1)];
                [m1, b1] = min(abs(dataBinsT{tr}-binWinT(ind,1)));
                [m2, b2] = min(abs(dataBinsT{tr}-binWinT(ind,2)));
                assert(all([m1 m2]<0.0001)); % bins should be almost identical
                binWin(ind,:) = [b1 b2];
                % average across bins within each window
                segments(:,tr,cati,ind) = cell2mat(cellfun(@(x) nanmean(x(binWin(ind,1):binWin(ind,2))), spk, 'UniformOutput', false))'; % ch x trial x cond x class (timebin)
            end
        end
        
        % average baselines across trials
        spk_baseline(:,cati) = squeeze(nanmean(chBaselines(:,:,cati),1));
        
        %one hot code the labels
        mod_1HC = zeros(size(segments,2)*(numBins), numBins);
        ct = 0; segTall = [];
        for tr = 1: size(segments,2) % for each trial
            for li = 1: (size(segments,4)) %for each label not baseline
                ct = ct+1;
                mod_1HC(ct, li)=1;
                segTall(:,ct) = segments(:,tr,cati,li);
            end
        end
        % append the ITI FR, duplicated to match the number of repetitions
        modBase_1HC = [mod_1HC; zeros(size(segments,2),numBins)];
        
        for ch = 1:128 %for each channel
            % append baseline FR
            segTallBase = [segTall(ch,:)'; ones(size(segments,2), 1) *spk_baseline(ch,cati)];
            mdl = fitlm(modBase_1HC,segTallBase);
            pLM_mod_byBin(ch,cati,:) = mdl.Coefficients.pValue(2:end); % pull p coefficients
        end
        % correct for multiple comparisons
        for te = 1:size(pLM_mod_byBin, 3) % for each variable
            pLM_mod_byBin(:, cati,te) = utils.bonf_holm(pLM_mod_byBin(:, cati,te));
        end
             
        % bootstrap
        disp(['beginning ' runNames{cati} ' tuning bootstrap'])
        for bt = 1:iter
            if rem(bt,200)==0
                disp(['on iteration ' num2str(bt) ' of ' num2str(iter)]);
            end
            
            % samples with replacement
            trToUseBoot = datasample(1:size(segments,2),size(segments,2)); 
            
            %one hot code the labels
            mod_1HC = zeros(size(segments,2)*(numBins), numBins);
            ct = 0; segTall = [];
            for tr = 1: size(segments,2) % for each trial
                for li = 1: (size(segments,4)) %for each label not baseline
                    ct = ct+1;
                    mod_1HC(ct, li)=1;
                    segTall(:,ct) = segments(:,trToUseBoot(tr),cati,li);
                end
            end
            
            % append the ITI FR, duplicated to match the number of repetitions
            modBase_1HC = [mod_1HC; zeros(size(segments,2),numBins)];
            for ch = 1:128 %for each channel
                % recompute baseline
                spk_baselineB = squeeze(nanmean(chBaselines(trToUseBoot,ch,cati),1));
                
                % append baseline FR
                segTallBase = [segTall(ch,:)'; ones(size(segments,2), 1) *spk_baselineB];
                mdl = fitlm(modBase_1HC,segTallBase);
                pLM_mod_byBin_boot(ch,cati,:,bt) = mdl.Coefficients.pValue(2:end); % not getting offset beta
            end
            % correct for multiple comparisons
            for te = 1:size(pLM_mod_byBin, 3) % for each variable
                pLM_mod_byBin_boot(:, cati,te,bt) = utils.bonf_holm(pLM_mod_byBin_boot(:, cati,te,bt));
            end
            
        end
    end
end

% compute CIs for the bootstrap
isSigBoot = squeeze(sum(pLM_mod_byBin_boot<0.05,1));
isTunedUpperCI = prctile(isSigBoot,97.5,3);
isTunedLowerCI = prctile(isSigBoot,2.5,3);

%% plots the number of tuned channels (Fig.5a)

alpha = 0.4; % transparency of plot
taskTimings = [0 0.5 1 1.5]; % in s, [start visual (down), start touch, end touch, end visual (up)]

% colormap
catHues = [0.364705882	0.733333333	0.278431373 % realistic
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
f1 = figure('Position',[200 500 950 375]);
for cati = 1:2
    subplot(1,2,cati)
    
    % extract values to plot where p-values are significant
    thisP = squeeze(pLM_mod_byBin(:, cati,:));
    sumP = sum(thisP<0.05);
    lowerCI = squeeze(isTunedLowerCI(cati,:));
    upperCI = squeeze(isTunedUpperCI(cati,:));
    upperP = upperCI-sumP;
    lowerP = sumP-lowerCI;
    
    % patches for when visual animations took place
    ylim([0 20]);yy = ylim();
    yp = [yy(1) yy(1) yy(2) yy(2)];
    patch([taskTimings(2:3) taskTimings(3:-1:2)]-0.5,yp,stimColor, 'FaceAlpha',0.5, 'LineStyle','none')
    hold on
    patch([taskTimings(1:2) taskTimings(2:-1:1)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    patch([taskTimings(3:4) taskTimings(4:-1:3)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    
    [h0, h1(cati)] = utils.plotsem_alpha_ul(binCenterTouch,sumP,...
        lowerP, upperP,catHues(cati,:),catHues(cati,:),alpha, '-', 3);
    set(h0,'edgec',catHues(cati,:));
    plot(binCenterTouch,sumP,'o', 'LineWidth',3,...
        'Color',catHues(cati,:), 'MarkerFaceColor',catHues(cati,:),'MarkerSize',4);

    title([runNames{cati}]);
    ylim([0 20])
    xlim([binCenterTouch(1) binCenterTouch(end)])
    xticks([-1:0.5:2])
end
sgtitle(['P1: S1 linear regression tuned channels relative to baseline'])

%% plots example channel activity (Fig.5b)
chs = [103 105]; % channels selected to plot
doSmooth = 1; % to perform smoothing on firing rates
smoothFactor = 5; % smoothing factor 
alpha = 0.4; % transparency of plotted lines

% colormap
catHues = [0.364705882	0.733333333	0.278431373 % realistic
    0.08627451	0.37254902	0.22745098]; % abstract

for chI = 1:numel(chs) % for each channel to plot

    ch = chs(chI);
    
    f2 = figure('Position', [150 90 480 425], 'Color',[1 1 1]);
    
    % patches for when visual animations took place
    ylim([0 4.5]);yy = ylim();
    yp = [yy(1) yy(1) yy(2) yy(2)];
    patch([taskTimings(2:3) taskTimings(3:-1:2)]-0.5,yp,stimColor, 'FaceAlpha',0.5, 'LineStyle','none')
    hold on
    patch([taskTimings(1:2) taskTimings(2:-1:1)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    patch([taskTimings(3:4) taskTimings(4:-1:3)]-0.5,yp,sideColor, 'FaceAlpha',0.5, 'LineStyle','none')
    
    % get firing rates and plot
    for cati = 2:-1:1
        data = (dataN(runCats{cati})); data = cat(1,data{:});
        dataAmps = cellfun(@(x) x.trialAmps, dataStruct(runCats{cati}), 'UniformOutput',0);
        dataAmps = cat(2,dataAmps{:});
        uVisStartT = cellfun(@(x) x.visStart, dataStruct(runCats{cati}), 'UniformOutput',0);
        uVisStartT = cat(2,uVisStartT{:});
        dataBinsT = cellfun(@(x) x.binsT, dataStruct(runCats{cati}), 'UniformOutput',0);
        dataBinsT = cat(2,dataBinsT{:});
        
        % extract catch trials only
        tmpInd = (dataAmps==0);
        data = data(tmpInd,:,:);
        uVisStartT = uVisStartT(tmpInd);
        dataBinsT = dataBinsT(tmpInd);
        
        % extract the bins for each trial
        binWin = []; spkSeg = [];
        for tr = 1: size(data,1)
            dataTr =  squeeze(data(tr,:));
            winTs = uVisStartT(tr)+[binTs(1) binTs(end)];
            
            [m1, b1] = min(abs(dataBinsT{tr}-winTs(1)));
            [m2, b2] = min(abs(dataBinsT{tr}-winTs(2)));
            assert(all([m1 m2]<0.0001)); % bins should be almost identical
            binWin(tr,:) = [b1 b2];
            spkSeg(tr,:,:) = cell2mat(cellfun(@(x) (x(binWin(tr,1):binWin(tr,2))'), dataTr, 'UniformOutput', false)); % tr x timebin x ch
        end
        spk = squeeze(spkSeg(:,:,ch));
         
        % smooth spike firing rates
        if doSmooth
            spk = sgolayfilt(spk,1,smoothFactor,[],2);
        end
        
        % extract data, removing nans
        a = 1:size(spk,2);
        b = nanmean(spk,1);
        c = 1.96*(nanstd(spk,[],1)/sqrt(size(spk,1)));      
        keepIndex = ~isnan(b) & ~isnan(c);
        a = a(keepIndex);
        b = b(keepIndex);
        c = c(keepIndex);
        
        % plot
        utils.plotsem_alpha_ul([binTs(1):0.05:binTs(end)]-0.5,b,c,c,catHues(cati,:),catHues(cati,:),alpha, '-',5);
    end
    yy = ylim;
    
    hold on
    xlim([binTs(1) binTs(end)]-0.5);
    title(['S1, ch ' num2str(ch)]);

    ylim([yy]);
end

%% plots the venn diagram of overlap of tuned channels between visual conditions (Fig.5c)

timesToTest = [-0.45:binWidth:0.95]; % all times when the visual movement was active

% get venn diagram values
tunedCh_r = []; tunedCh_a = [];
for ti = 1:numel(timesToTest)
    [~,binToTest] = min(abs(binCenterTouch-timesToTest(ti)));
    tunedCh_r = [tunedCh_r; find(squeeze(pLM_mod_byBin(:, 1,binToTest))<0.05)];
    tunedCh_a = [tunedCh_a; find(squeeze(pLM_mod_byBin(:, 2,binToTest))<0.05)];
end
tunedCh_r = unique(tunedCh_r);
tunedCh_a = unique(tunedCh_a);

disp('P1: tuning results')
disp(['realistic: number of channels tuned in visual bins: ' num2str(numel(tunedCh_r))])
disp(['abstract: number channels tuned in visual bins: ' num2str(numel(tunedCh_a))])

% plot venn diagram
numBoth = numel(intersect(tunedCh_r,tunedCh_a));
numRonly = numel(tunedCh_r)-numBoth;
numAonly = numel(tunedCh_a)-numBoth;
numBoth + numRonly + numAonly;
figure;
utils.vennX([numRonly numBoth numAonly],0.005); % this is just for visualization, make pretty in affinity
v1 = gcf;
title('P1: Realistic vs abstract venn diagram of tuned channels')
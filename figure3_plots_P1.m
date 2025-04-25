% figure3_plots_P1

% for participant P1

% plots the intensity ratings of ICMS sensations (Fig.3a)
% plots the locations of ICMS sensations (Fig.3b)
% plots descriptor words of ICMS sensations(Fig.3c)

% Isabelle Rosenthal 2025

% load data
load(['..\Data\stimContext_Timing_P1_preprocessedSpks.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots the intensity ratings of ICMS sensations (Fig.3a)

ampsUsed = [100 60 30 0]; % current amplitudes tested
timesUsed = [-300  -150     0   150   300 Inf]; % timing offsets tested for P1
runNames = {'realistic','abstract'};

ampMap = [4 3 2 1]; % order the amplitudes the way we want

% colormap
ampColors = [57 57 57; % 100
143 143 143; % 60
227 227 227]/255; % 100

fbt = 0; % felt sensations in baseline counter
ft = 0; % felt sensations visual conditions counter
for di = 1:size(dataStruct,1) % for each day
    for si = 1:size(dataStruct,2) % for each run
        runType = dataStruct{di,si}.runType; % 1=realistic, 2=abstract
        
        % get the trials where a sensation was elicited
        feltTrials = ~isnan(dataStruct{di,si}.trialIntensity);
        feltInds = find(feltTrials);
        
        if runType>0 % if not baseline
            for fi = 1:numel(feltInds) % get the felt intensities
                ft = ft+1;
                intVec(ft) = dataStruct{di,si}.trialIntensity(feltInds(fi));
                condVec{ft} = runNames{runType};
                ampVec{ft} = num2str(dataStruct{di,si}.trialAmps(feltInds(fi)));
                timeVec(ft) = dataStruct{di,si}.trialVisRelativeToStim(feltInds(fi));
            end
        else
            for fi = 1:numel(feltInds) % get baseline felt intensities too
                fbt = fbt+1;
                intBaseVec(fbt) = dataStruct{di,si}.trialIntensity(feltInds(fi));
                ampBaseVec{fbt} = num2str(dataStruct{di,si}.trialAmps(feltInds(fi)));
            end
        end
    end
end

disp('P1 results: sensation intensities')
% 2-way anova on condition x currrent amp, within baseline, realistic, abstract
allInt = [intBaseVec intVec];
allCond = cell(1,size(allInt,2));
allCond(1:size(intBaseVec,2)) = {'baseline'};
allCond(size(intBaseVec,2)+1:end) = condVec;
allAmp = cell(1,size(allInt,2));
allAmp(1:size(intBaseVec,2)) = ampBaseVec;
allAmp(size(intBaseVec,2)+1:end) = ampVec;
p = anovan(allInt,{allCond, allAmp},'model','interaction',...
    'varnames',{'condition','current'});
disp('2-way anova on intensity ratings')
disp(['effect of condition: p=' num2str(p(1))]);
disp(['effect of current: p=' num2str(p(2))]);
disp(['effect of condition*current: p=' num2str(p(3))]);

disp('*************************')

% 3 way anova
% condition x current amp x timing offset within realistic vs abstract trials?
p = anovan(intVec,{condVec, ampVec,timeVec},'model','interaction',...
    'varnames',{'condition','current','time'});
disp('3-way anova on intensity ratings within realistic vs abstract trials')
disp(['effect of condition: p=' num2str(p(1))]);
disp(['effect of current: p=' num2str(p(2))]);
disp(['effect of timing offset: p=' num2str(p(3))]);
disp(['effect of condition*current: p=' num2str(p(4))]);
disp(['effect of condition*time: p=' num2str(p(5))]);
disp(['effect of current*time: p=' num2str(p(6))]);

disp('*************************')

% unpaired ttest: intensities between 100 vs 60uA
[h p] = ttest2(intVec(strcmp(ampVec,'100')),intVec(strcmp(ampVec,'60')));
disp(['unpaired ttest of intensities between 100 and 60: p=' num2str(p)]);

% do a violin plot
figure('Position',[900 900 1500 350]);
utils.violinplot(intVec, ampVec, 'GroupOrder', {'100','60','30'},...
    'HalfViolin', 'left','QuartileStyle','none','ViolinAlpha',1,'DataStyle','histogram',...
    'ViolinColor',ampColors);
title(['P1 sensation intensity reports (realistic and abstract trials)'])
ylabel('sensation intensity')
ylim([0.5 8.5])
xlabel('ICMS current amplitude')

% do a violin plot for the baseline
figure('Position',[900 900 1500 350]);
utils.violinplot(intBaseVec, ampBaseVec, 'GroupOrder', {'100','60','30'},...
    'HalfViolin', 'left','QuartileStyle','none','ViolinAlpha',1,'DataStyle','histogram',...
    'ViolinColor',ampColors);
title(['P1: sensation intensity reports (baseline trials)']);
ylabel('sensation intensity')
ylim([0.5 8.5])
xlabel('ICMS current amplitude')

%% plots the locations of ICMS sensations (Fig.3b)
% generates one plot for each condition and current amplitude tested

% background image
bodyImg = 'imgs\stimLocationArms.png';

ampsUsed = [100 60 30 0]; % current amplitudes tested
runNames = {'baseline','realistic','abstract'};

% locations reported and their grid coordinates
locs =[ {'c8'}    {'c9'}   {'d9'}  {'d12'}    {'d13'}    {'e13'}   ...
    {'w8'}  {'w9'}  {'x10'}    {'x8'}    {'x9'}    {'x7'}];
coords = [     0.1150   -0.7400 %c8
    0.1150   -0.708 %c9
    0.1460   -0.708 %d9
    0.1460   -0.615 %d12
    0.1460   -0.5820 %d13
    0.1750   -0.5820 %e13
    0.886   -0.74 %w8
    0.886   -0.708 %w9
    0.918   -0.6755 %x10
    0.918   -0.7400 %x8
    0.918   -0.708 %x9
    0.918   -0.7700]; %x7

% colormap
cmap = (hot(125));
cmap = cmap(1:100,:);

coordCountByCA = zeros(numel(locs),3, 3); % count locations by condition and amplitude
for di =  1:size(dataStruct,1) % for each day
    for si = 1:size(dataStruct,2) % for each run
        runType = dataStruct{di,si}.runType;
        feltTrials = ~isnan(dataStruct{di,si}.trialIntensity);
        feltInds = find(feltTrials); % was there a sensation elicited
        trialAmps = dataStruct{di,si}.trialAmps(feltInds); % current amplitude used
        
        % count the locations reported
        for trI = 1:numel(feltInds)
            cI = strcmp(dataStruct{di,si}.trialLoc(feltInds(trI)), locs);
            ampU = find(ampsUsed==trialAmps(trI));
            coordCountByCA(cI, runType+1,ampU) =  coordCountByCA(cI, runType+1,ampU) + 1;      
        end
        
    end
end

% plot the locations of ICMS sensations by current and condition
totalCountByCA = sum(coordCountByCA,1); % sum within cond and amp
perCountCA = 100*coordCountByCA./totalCountByCA;
for rI = 1:size(coordCountByCA,2) % for every run type
    for aI = 1:size(coordCountByCA,3) % for every amplitude tested
        figure('Position',[100 100 880 775], 'Color', [1 1 1]);
        axis tight;
        hold on
        I = imread(bodyImg);
        h = image(xlim,-ylim,I);
        uistack(h,'bottom')
        for lI = 1:numel(locs)
            if coordCountByCA(lI,rI,aI)>0
                plot(coords(lI,1),coords(lI,2),'s','MarkerEdgeColor','k',...
                    'MarkerFaceColor',cmap(round(perCountCA(lI,rI,aI)),:), 'MarkerSize',19);
            end
        end
        colormap(cmap);
        colorbar;
        title(['P1: percentage of sensations in ' runNames{rI} ' trials, amp=' num2str(ampsUsed(aI))])
    end
end

%% plots descriptor words of ICMS sensations(Fig.3c)

timesUsed = [-300  -150     0   150   300]; % timing offsets tested for P1
ampsUsed = [100 60 30]; % current amplitudes tested
runNames = {'baseline','abstract','realistic'};
runMap = [1 3 2]; % remap conditions how we want: baseline abstract realistic

% unique words used to describe ICMS sensations (participant generated)
uniqueWds = {'touch','pinch','squeeze','grab','move right','move left','move up'};

% colormap
wordHues = [0.803921569	0.603921569	0.968627451
1	0.690196078	0
0.996078431	0.380392157	0
0.862745098	0.149019608	0.498039216
0.392156863	0.560784314	1
0.223529412	0.333333333	0.619607843
0.431372549	0.701960784	0.552941176];

% count the words used
wordsDesc = {}; % list of words used by session and condition
wordsDescA = {}; % list of words used by session, condition, and amplitude
wdCt_byA = []; % count of words used by session, condition, and amplitude
wdCt_justA = []; % count of words used by condition and amplitude
for di = 1:size(dataStruct,1) % for each day
    for si = 1:size(dataStruct,2) % for each run
        runType = dataStruct{di,si}.runType;
        feltTrials = ~isnan(dataStruct{di,si}.trialIntensity);        
        feltInds = find(feltTrials);
        
        % get the word descriptors used
        wordsDesc{di,runMap(runType+1)} = [dataStruct{di,si}.trialDescription(feltInds)];
        for ai = 1:3
            thisA = find((dataStruct{di,si}.trialAmps==ampsUsed(ai))& feltTrials);
            wordsDescA{di,runMap(runType+1),ai} = [dataStruct{di,si}.trialDescription(thisA)];
        end
    end
    
    % get the word count by amplitude
    for rt = 1:3 % for each condition
        for uw = 1:numel(uniqueWds) % for each word
            for ai = 1:3 % for each amplitude
                wdCt_byA(di,rt,ai,uw) = mean(strcmp(wordsDescA{di,rt,ai}, uniqueWds{uw}))*100; % this is percentage of words per session in that amp
            end
        end
    end
    % get the word count combined by days - for sign rank test across conditions
    wdCt_justA(di,:,:) = squeeze(sum(wdCt_byA(di,:,:,:),2));
end
wdCt_allA = squeeze(nanmean(wdCt_byA,1)); % take the mean of the percentages across sessions

% check for words not in list
newUnique = unique(horzcat(wordsDesc{:}));
if numel(newUnique)~=numel(uniqueWds)
    error('WORDS NOT IN LIST')
end

% plot descriptor words of ICMS sensations(Fig.3c)
% make pie charts for every condition and amplitude
figure('Position',[200 300 1150 900]);
for runType=1:3 % for every condition
    for ai = 1:2 % for every current amplitude
        subplot(3,2,((runType-1)*2) + ai)
        % only plot words that had at least one occurance
        valsToPlot =  squeeze(wdCt_allA(runType,ai,:));
        indsToPlot = valsToPlot~=0;
        valsToPlot = valsToPlot(indsToPlot);
        
        p = pie( valsToPlot);
        ax = gca(); ax.Colormap = wordHues(indsToPlot,:);
        
        % get number of words used overall for the title
        numWordsA = sum(cellfun(@(x) size(x,2), wordsDescA(:,runType,ai)));
        
        title([runNames{runType} ': ' num2str(ampsUsed(ai))  'uA (n=' num2str(numWordsA) ')']);
        
        legend(uniqueWds(indsToPlot), 'Location','westoutside')     
    end
end
sgtitle('P1: ICMS-elicited sensation descriptors')

disp('P1 results: sensation descriptors')
[pTouch_100(1)] = signrank(wdCt_byA(:,1,1,1),wdCt_byA(:,3,1,1)); % baseline vs real touch occurrence
[pTouch_100(2)] = signrank(wdCt_byA(:,1,1,1),wdCt_byA(:,2,1,1)); % baseline vs abstract touch occurrence
[pTouch_100(3)] = signrank(wdCt_byA(:,3,1,1),wdCt_byA(:,2,1,1)); % baseline vs abstract touch occurrence
pTouch_100 = utils.bonf_holm(pTouch_100); % bonferroni correct  
disp(['sign rank 100A baseline vs real ''touch'' occurence: ' num2str(pTouch_100(1))])
disp(['sign rank 100A baseline vs abstract ''touch'' occurence: ' num2str(pTouch_100(2))])
disp(['sign rank 100A realistic vs abstract ''touch'' occurence: ' num2str(pTouch_100(3))])
[pTouch_60(1)] = signrank(wdCt_byA(:,1,2,1),wdCt_byA(:,3,2,1)); % baseline vs real touch occurrence
[pTouch_60(2)] = signrank(wdCt_byA(:,1,2,1),wdCt_byA(:,2,2,1)); % baseline vs abstract touch occurrence
[pTouch_60(3)] = signrank(wdCt_byA(:,3,2,1),wdCt_byA(:,2,2,1)); % baseline vs abstract touch occurrence
pTouch_60 = utils.bonf_holm(pTouch_60);% bonferroni correct
disp(['sign rank 60A baseline vs real ''touch'' occurence: ' num2str(pTouch_60(1))])
disp(['sign rank 60A baseline vs abstract ''touch'' occurence: ' num2str(pTouch_60(2))])
disp(['sign rank 60A realistic vs abstract ''touch'' occurence: ' num2str(pTouch_60(3))])
  
% does touch get used more at 60 than at 100? Looking across conditions
[pTouch_60v100_allCond] = signrank(wdCt_justA(:,1,1),wdCt_justA(:,2,1));  % 100a touch vs 60a touch
disp(['sign rank 60A vs 100A across conditions ''touch'' occurence: ' num2str(pTouch_60v100_allCond)])

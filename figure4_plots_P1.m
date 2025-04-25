% figure4_plots_P1

% for participant P1

% plots the temporal binding window behavioral data for 100uA and 60uA
% trials (Fig.4a, 4c)
% plots the Gaussian fit and parametric bootstrap of the temporal binding
% window for 100uA and 60uA trials (Fig.4b, 4c)

% Isabelle Rosenthal 2025

% if you want to load in the bootstrap analysis with 5000 iterations
% instead of rerunning it yourself
toLoad = 1; % 1=load in prerun analysis, 0=run it yourself
% if you run yourself, bootstrapped values may be slightly different due to
% change in the random seed

% load data
load(['..\Data\stimContext_Timing_P1_preprocessedSpks.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots the temporal binding window behavioral data
% for 100uA and 60uA trials (Fig.4a, 4c)

ampsUsed = [100 60 30 0]; % current amplitudes tested
timesUsed = [-300  -150     0   150   300 Inf]; % timing offsets tested for P1
runNames = {'realistic','abstract'};

% colormap
timingReportColors = [...
    [255 127 101];...   % vision first
    [0 0 0];...   % simultaneous
    [86 158 255];...   % sensation first
    ]/255;

offsetMat = zeros(5,3,2,4); %count of reported order over all sessions, times x ans x runtype x amplitude
totalAnsByTime = zeros(9,5,2); % total felt sensations, days x times x runType
StimFirstAnsByTimeByAmp = zeros(9,5,4,2); % count reported order by session, days x times x amp x runtype
VisFirstAnsByTimeByAmp = zeros(9,5,4,2); % count reported order by session, days x times x amp x runtype
sameAnsByTimeByAmp = zeros(9,5,4,2); % count reported order by session, days x times x amp x runtype
totalAnsByTimeByAmp =  zeros(9,5,4,2); % days x times x amp x runtype
ft = 0; % felt trial counter
isSync = []; c_amps = []; c_time = []; c_cond = []; % aggregate data together across all runs
for di = 1:size(dataStruct,1) % for each day
    for si = 1:size(dataStruct,2) % for each run
        if dataStruct{di,si}.runType ~=0
            runType = dataStruct{di,si}.runType;
            feltTrials = ~isnan(dataStruct{di,si}.trialIntensity);
            feltInds = find(feltTrials);
            senAns = dataStruct{di,si}.trialSenOrder;
            senOrder = nan(size(senAns));
            senOrder(strcmp(senAns,'Vision First')) = -1;
            senOrder(strcmp(senAns,'Stim First')) = 1;
            senOrder(strcmp(senAns,'Simultaneous')) = 0;
            
            timesUsed = unique(dataStruct{di,si}.trialVisRelativeToStim); % inf is catch trial
            
            for trI = 1:numel(feltInds)
                aI = dataStruct{di,si}.trialAmps(feltInds(trI))==ampsUsed;
                sI = dataStruct{di,si}.trialVisRelativeToStim(feltInds(trI))==timesUsed;
                if strcmp(dataStruct{di,si}.trialSenOrder(feltInds(trI)),'Vision First')
                    offsetMat(find(sI),1,runType, aI) =  offsetMat(find(sI),1,runType, aI)+1;
                    VisFirstAnsByTimeByAmp(di,sI, aI, runType) = VisFirstAnsByTimeByAmp(di,sI, aI, runType) + 1;
                elseif strcmp(dataStruct{di,si}.trialSenOrder(feltInds(trI)),'Simultaneous')
                    offsetMat(find(sI),2,runType, aI) =  offsetMat(find(sI),2,runType, aI)+1;
                    sameAnsByTimeByAmp(di,sI, aI, runType) = sameAnsByTimeByAmp(di,sI, aI, runType) + 1;
                elseif strcmp(dataStruct{di,si}.trialSenOrder(feltInds(trI)),'Stim First')
                    offsetMat(find(sI),3,runType, aI) =  offsetMat(find(sI),3,runType, aI)+1;
                    StimFirstAnsByTimeByAmp(di,sI, aI, runType) = StimFirstAnsByTimeByAmp(di,sI, aI, runType) + 1;
                end
                
                totalAnsByTime(di, sI, runType) = totalAnsByTime(di, sI, runType) + 1;
                totalAnsByTimeByAmp(di, sI,  aI, runType) = totalAnsByTimeByAmp(di, sI,  aI, runType) + 1;
                
                % get data for logistic regression test
                ft = ft+1;
                isSync(ft) = (senOrder(feltInds(trI))==0)+1;
                c_amps(ft) =  dataStruct{di,si}.trialAmps(feltInds(trI));
                c_time(ft) =  dataStruct{di,si}.trialVisRelativeToStim(feltInds(trI));
                c_cond(ft) =  dataStruct{di,si}.runType;
            end
        end
    end
end

disp('P1 results: temporal order judgements')

perReportVisFirst = []; perReportStimFirst = []; perReportSame = []; % percentage this answer was reported
perVisF = []; perSame = []; perStimF = []; % percentage reported within sessions
perVisFSem = []; perSameSem = []; perStimFSem = []; % std of percentage across sessions
for aI = 1:2 % for 100uA, 60ua
    perReportVisFirst(:,:,aI) = squeeze(100*offsetMat(:,1,:,aI)./sum(offsetMat(:,:,:,aI),2)); % sum across amplitudes
    perReportStimFirst(:,:,aI) =  squeeze(100*offsetMat(:,3,:,aI)./sum(offsetMat(:,:,:,aI),2));
    perReportSame(:,:,aI) =  squeeze(100*offsetMat(:,2,:,aI)./sum(offsetMat(:,:,:,aI),2));
    
    perVisF(:,:,:,aI) = squeeze(100*VisFirstAnsByTimeByAmp(:,:,aI,:)./totalAnsByTimeByAmp(:,:,aI,:));% days x times x runType
    perVisFSem(:,:,aI) = squeeze(nanstd(perVisF(:,:,:,aI))/sqrt(size(perVisF(:,:,:,aI),1)));
    perSame(:,:,:,aI) = squeeze(100*sameAnsByTimeByAmp(:,:,aI,:)./totalAnsByTimeByAmp(:,:,aI,:));% days x times x runType
    perSameSem(:,:,aI) = squeeze(nanstd(perSame(:,:,:,aI))/sqrt(size(perSame(:,:,:,aI),1)));
    perStimF(:,:,:,aI) = squeeze(100*StimFirstAnsByTimeByAmp(:,:,aI,:)./totalAnsByTimeByAmp(:,:,aI,:));% days x times x runType
    perStimFSem(:,:,aI) = squeeze(nanstd(perStimF(:,:,:,aI))/sqrt(size(perStimF(:,:,:,aI),1)));
    
    % test if the area under the curve is the same between abstract and realistic trials
    perSame_nn = perSame(:,:,:,aI);
    perSame_nn(isnan(perSame(:,:,:,aI)))=0;
    areaReal = trapz(perSame_nn(:,:,1),2);
    areaAbs = trapz(perSame_nn(:,:,2),2);
    [~, p] = ttest(areaReal,areaAbs, 'Tail','right');
    disp([num2str(ampsUsed(aI)) 'uA: is there is a diff between realistic and abstract area under "same" curve?'])
    disp(['one-tailed paired ttest p= ' num2str(p)])
    
    % plot the same plots but separated by amplitude
    f1 = figure('Position',[400 350 1000 400]);
    for rI = 1:2
        subplot(1,2,rI)
        plot([0 0],[0 100],'k:')
        hold on
        % plot 25% line
        plot([-320 320],[25 25],'k:')
        a(1) = plot(timesUsed(1:5),perReportVisFirst(:,rI,aI),'LineWidth', 3, 'Color',timingReportColors(1,:));
        hold on
        a(2) = plot(timesUsed(1:5),perReportSame(:,rI,aI),'LineWidth', 3,'Color',timingReportColors(2,:));
        a(3) = plot(timesUsed(1:5),perReportStimFirst(:,rI,aI),'LineWidth', 3,'Color',timingReportColors(3,:));
        errorbar(timesUsed(1:5),perReportVisFirst(:,rI,aI), perVisFSem(:,1,aI),perVisFSem(:,rI,aI),'.', 'Color',timingReportColors(1,:),...
            'LineWidth',2);
        errorbar(timesUsed(1:5),perReportSame(:,rI,aI), perSameSem(:,1,aI),perSameSem(:,rI,aI),'.', 'Color',timingReportColors(2,:),...
            'LineWidth',2);
        errorbar(timesUsed(1:5),perReportStimFirst(:,rI,aI), perStimFSem(:,1,aI),perStimFSem(:,rI,aI),'.', 'Color',timingReportColors(3,:),...
            'LineWidth',2);
        legend(a, 'Vis First','Same','Stim First', 'Location','northwest')
        xlabel('vis relative to stim actual times')
        ylabel('% reported')
        title(runNames{rI})
        ylim([0 100])
        xlim([-320 320])
    end
    sgtitle(['P1: ' num2str(ampsUsed(aI)) 'uA trials'])
end

disp('**********************************')
% logistic regression test
c_X = [c_amps; c_time; c_cond]';
[B_d, dev_d, stats_d] = mnrfit(c_X, isSync', 'model','nominal');
p_d = [stats_d.p(2:4)];
% bonferroni correct since we are combining/comparing categories
p_d = utils.bonf_holm(p_d);
% the first term in B is the intercept
disp(['logistic regression test'])
disp(['significance of amp in giving ''simultaneous'' answer: ' num2str(p_d(1))])
disp(['significance of timing in giving ''simultaneous'' answer: ' num2str(p_d(2))])
disp(['significance of condition in giving ''simultaneous'' answer: ' num2str(p_d(3))])

%% run the Gaussian fit and parametric bootstrap analysis
% of the temporal binding window for 100uA and 60uA trials
% OR load in the prerun file

iter = 5000; % number of bootstrap iterations to run
xVals = -550:0.1:550; % the offset times to fit to
catNames = {'realistic','abstract'}; % condition names
ampNames = {'100','60'}; % current amplitude names

disp('**********************************')

if toLoad == 1
    % if you're loading in former results
    load(['..\Data\stimContext_Timing_P1_gaussianBootByAmp.mat']);
    disp('results loaded successfully')
    disp('**********************************')
else
    
    gausParamsA = []; % Gaussian parameters
    gausValsA = []; % % values from the Gaussian fit
    gausParamsA_boot = []; % Gaussian parameters, bootstrapped
    gausValsA_boot = []; % values from the Gaussian fit, bootstrapped
    
    sameAnsPerByAmp = sameAnsByTimeByAmp./totalAnsByTimeByAmp; % totalAnsByTime = total number of trials per session
    % get rid of nans
    sameAnsPerByAmp(isnan(sameAnsPerByAmp))=0;
    
    disp('Note that bootstrapped values may vary slightly due to randomness in the distribution.')
    disp('Set a random seed to ensure replicability.')
    
    for cati = 1:2 % for each condition
        for ai = 1:2 % for 100uA and 60uA
            
            thisPer = sameAnsPerByAmp(:,:,ai,cati); % percent time simultaneous was the answer, divided by session
            thisPer = thisPer(:); % collapse across sessions and times tested
            thisX = repmat(timesUsed(1:5),size(sameAnsPerByAmp,1),1); % make matrix of times tested
            thisX = thisX(:); % collapse across sessions and times tested
            
            % fit a gaussian to all data: fit to the probability on any given sess
            % that sync was the answer
            options = fitoptions('gauss1', 'StartPoint',[0.5 0 50],'Lower', [0 -Inf -Inf], 'Upper', [1 Inf Inf]);
            [gaufit, gofit]= fit(thisX, thisPer,'gauss1', options); % fit one gaussian to the curve
            gausParamsA(cati,ai,:) = [ gaufit.a1  gaufit.b1  gaufit.c1]; % peak height, peak time, gausSTD
            
            % get fitted line from formula which is 'formula(gaufit)'
            gausValsA(cati,ai,:) = gausParamsA(cati,ai,1).*exp(-((xVals'-gausParamsA(cati,ai,2))./gausParamsA(cati,ai,3)).^2);
            
            % get the values specifically at the times tested, for the bootstrap
            pOffsetA =  gausParamsA(cati,ai,1).*exp(-((timesUsed(1:5)'-gausParamsA(cati,ai,2))./gausParamsA(cati,ai,3)).^2);
            
            % bootstrap
            % iterate 'iter' times over the same number of trials from the fitted
            % gaussian, using a binomial distribution, and fit a gaussian to each
            % one of these randomly generated sets of data
            disp(['beginning ' catNames{cati} ', ' ampNames{ai} 'A parametric bootstrap'])
            for bi = 1:iter
                if mod(bi,500)==0
                    disp(['on boot ' num2str(bi) '/' num2str(iter)])
                end
                
                % generate new trials in the same quantity as the original data
                sameAnsPerA_boot = [];
                for sesh = 1:size(sameAnsByTimeByAmp,1) % for each session
                    for ti = 1:size(sameAnsByTimeByAmp,2) % for each offset
                        sameAnsByTimeA_boot = binornd(totalAnsByTimeByAmp(sesh,ti,ai,cati),pOffsetA(ti)); %n=total trials, p=from gaussian fit to original data
                        sameAnsPerA_boot(sesh,ti) = sameAnsByTimeA_boot/totalAnsByTimeByAmp(sesh,ti,ai,cati);
                    end
                end
                % get rid of nans
                sameAnsPerA_boot(isnan(sameAnsPerA_boot))=0;
                
                % fit the gaussian to this iteration
                sameAnsPerA_boot = sameAnsPerA_boot(:);
                [gaufit, gofit]= fit(thisX, sameAnsPerA_boot,'gauss1', options); % fit one gaussian to the curve
                gausParamsA_boot(cati,ai,:,bi) = [ gaufit.a1  gaufit.b1  gaufit.c1]; % peak height, peak time, gausSTD
                gausValsA_boot(cati,ai,:,bi) = gausParamsA_boot(cati,ai,1,bi).*exp(-((xVals'-gausParamsA_boot(cati,ai,2,bi))./gausParamsA_boot(cati,ai,3,bi)).^2); % get fits from formula which is formula(gausfit{cati})
            end
            
            % get 95% CIs on the gaussian fit as well as on peak, time, and std
            gausValsA_CI05(cati,ai,:,:) = [prctile(squeeze(gausValsA_boot(cati,ai,:,:))',2.5); prctile(squeeze(gausValsA_boot(cati,ai,:,:))',97.5)];
            gausParamsA_CI05(cati,ai,:,:) = [prctile(squeeze(gausParamsA_boot(cati,ai,:,:))',2.5); prctile(squeeze(gausParamsA_boot(cati,ai,:,:))',97.5)];
        end
    end
end

%% plots the Gaussian fit and parametric bootstrap of the temporal binding
% window for 100uA and 60uA trials (Fig.4b, 4c)

alpha = 0.3; % transparency level in the plots

% colormap
condsHues = [93	187	71; % realistic
    22	95	58]/255; % abstract

% initialize variables for calulating the 25% points on the Gaussian
x25sABoot = [];
xVals_added = -550:0.1:800; % extended x vals for plotting

disp('**********************************')

for ai = 1:2 % for 100, 60uA trials
    figure('Position',[700 500 930 500]);
    hold on
    for cati = 1:2 % for each condition
        % plot the bootstrapped distribution
        [~, h1(cati)] = utils.plotsem_alpha_ul(xVals,squeeze(gausValsA(cati,ai,:))',...
            (squeeze(gausValsA(cati,ai,:))-squeeze(gausValsA_CI05(cati,ai,1,:)))',...
            (squeeze(gausValsA_CI05(cati,ai,2,:))-squeeze(gausValsA(cati,ai,:)))',condsHues(cati,:),condsHues(cati,:),alpha, '-', 2);
        
        % plot PSS
        errorbar(gausParamsA(cati,ai,2), gausParamsA(cati,ai,1),gausParamsA(cati,ai,2)-gausParamsA_CI05(cati,ai, 1,2),...
            gausParamsA_CI05(cati,ai, 2,2)-gausParamsA(cati,ai,2), 'horizontal','ko','MarkerFaceColor','k')
        plot(gausParamsA(cati,ai,2), gausParamsA(cati,ai,1),'ko','MarkerFaceColor','k')
        
        % calculate the 25% points on either side of the Gaussian
        [~,peakInd] = min(abs(gausValsA(cati,ai,:)-gausParamsA(cati,ai,1)));
        [~,L25] = min(abs(gausValsA(cati,ai,1:peakInd)-0.25));
        [~,R25] = min(abs(gausValsA(cati,ai,peakInd+1:end)-0.25));
        R25 = R25+peakInd;
        x25sA(cati,ai,:) = xVals([L25 R25]);
 
        for iter = 1:size(gausValsA_boot,4)
            [~,peakIndB] = min(abs(gausValsA_boot(cati,ai,:,iter)-gausParamsA_boot(cati,ai,1,iter)));
            [~,L25B] = min(abs(gausValsA_boot(cati,ai,1:peakInd,iter)-0.25));
            [~,R25B] = min(abs(gausValsA_boot(cati,ai,peakInd+1:end,iter)-0.25));
            R25B = R25B+peakIndB;
            x25sABoot(cati,ai,:,iter) = xVals_added([L25B R25B]);
        end
    end
    xlim([xVals(1) xVals(end)])
    ylim([0 1])
    % plot JND line
    plot([xlim],[0.25 0.25],'k:')
    % plot 0 time line
    plot([0 0],[ylim],'k:')
    legend(h1,'realistic','abstract', 'Location','northwest')
    title(['P1: ' ampNames{ai} 'uA, Gaussian fit to simultaneous answer rates'])
    
    % get PSS (peak of each gaussian)
    disp(['P1: ' ampNames{ai} 'uA Gaussian values'])
    disp('PSS (point of subjective simultaneity/peak of each gaussian)')
    disp(['realistic PSS = ' num2str(gausParamsA(1,ai, 2)) 'ms, 95%CI [' num2str(gausParamsA_CI05(1,ai, :,2)) ']'])
    disp(['abstract PSS = ' num2str(gausParamsA(2, ai,2)) 'ms, 95%CI [' num2str(gausParamsA_CI05(2, ai,:,2)) ']'])
    
    for cati = 1:2
        bSidedA(1,:) = (gausParamsA_boot(cati,ai,2,:)-x25sABoot(cati,ai,1,:));
        bSidedA(2,:) = (gausParamsA_boot(cati,ai,2,:)-x25sABoot(cati,ai,2,:));
        JNDA_boot(cati,ai,:) = mean(abs(bSidedA));
    end
    JNDA_CIs(ai,:,:) = [prctile(JNDA_boot(:,ai,:), 2.5,3) prctile(JNDA_boot(:,ai,:), 97.5,3)]; % ai x cati x prc
    
    % get average JND
    disp('JND (just noticeable difference) aka 0.25 point of gaussians')
    disp(['realistic avg JND: ' num2str(mean(abs(gausParamsA(1,ai, 2)-x25sA(1,ai,:)))),...
        ' [' num2str(JNDA_CIs(ai,1,:)) ']' ])
    disp(['abstract avg JND: ' num2str(mean(abs(gausParamsA(2,ai, 2)-x25sA(2,ai,:)))),...
        ' [' num2str(JNDA_CIs(ai,2,:)) ']' ])
end

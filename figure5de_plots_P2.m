% figure5de_plots_P2

% for participant P2

% plots the RDMs of neural activity (Fig. 5d)
% plots the MDS diagrams of neural activity (Fig. 5e)

% Isabelle Rosenthal 2025

% load in the save file of RDM analysis
% the save file is already provided, but can be regenerated via the 
% python script figure5de_constructRDMs_P2.py
pyRDM = load(['..\Data\stimContext_Timing_P2_RDMs.mat']);
% if you run yourself, RDM values may be slightly different due to
% change in the random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots the RDMs of neural activity (Fig. 5d)

condRDM = pyRDM.RDM; % RDM
condNames = {'ITI_R', 'ITI_A', 'Down_R', 'Down_A', 'Touch_R', 'Touch_A',...
    'Up_R', 'Up_A'}; % phases of the task
condIDs = [1 2 3 4 5 6 7 8];
assert(all(condIDs==pyRDM.condLabelOrder.condLabels))

% visualize the RDM
close all % MDS code can mess with figures so close all figures first
numConds = numel(condIDs);

figure('Position',[50 200 1000 650]);
cmap = hot(100);

imagesc(condRDM);
yticks(1:1:numConds);
xticks(1:1:numConds);
xticklabels(condNames);
yticklabels(condNames);
xtickangle(-90);

caxis([0 0.06]);
colormap(cmap);
cc = colorbar('northoutside');
axis square
axis tight

title(['P2: S1 timing sorted RDM, CV noise-corr Mahalanobis']);

% are neural patterns more similar by phase than by conditions?
withinPAcrossCDist = [condRDM(1,2); condRDM(3,4); condRDM(5,6);condRDM(7,8);]; %ITI/ITI, Down/Down, etc
withinCAcrossPDist = [condRDM(1,3); condRDM(1,5);condRDM(1,7);... %ITIr/Downr, ITIr/Downr, ITIr/Upr within cond/across phase
    condRDM(2,4);condRDM(2,6);condRDM(2,8);... %ITIa/Downa, ITIa/Toucha, ITIa/Upa
    condRDM(3,5);condRDM(3,7);...
    condRDM(4,6);condRDM(4,8);...
    condRDM(5,7);...
    condRDM(6,8);];
[p h]=ranksum(withinPAcrossCDist,withinCAcrossPDist); %unpaired ttest
disp(['P2: rank sum test on distances within phases/across conds vs distances within conds/across phase:'])
disp(['p=' num2str(p)])

%% plots the MDS diagrams of neural activity (Fig. 5e)
% using toolbox from https://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/

close all % MDS code can mess with figures so close all figures first

userOptions = rsa.rsa_defineUserOptions_iar('StimContext_Timing');

% establish plotting conventions
localOptions.fontSize = 10;
userOptions.dotSize=20;
rdmStruct.name = 'StimContext_Timing';
localOptions.figureNumber = 1;
rdmStruct.color = [ 0 0 0];
phaseMarker = {'s','v','o','^',};
phaseMap = [1 1 2 2 3 3 4 4];
catMap = [1 2 1 2 1 2 1 2];

% colormap
condsHues = [93	187	71; % realistic
    22	95	58]/255; % abstract

userOptions.conditionColours =[];userOptions.markerStyles =[];
for mod = 1:numel(condIDs)
    ph = phaseMap(condIDs(mod)); cat = catMap(condIDs(mod));
    userOptions.conditionColours(mod,:)=condsHues(cat,:);
    userOptions.markerStyles{mod} = phaseMarker{ph};
end

localOptions.rotateAngle = 30*pi/180; % angle to rotate each bin, customized. Should be with touches pointing to northeast

cf = condRDM;
% edit to not plot any negative values
numNeg = (sum(sum(cf<0))-8)/2; % divide because symmetric and remove diagonal
cf(cf<=0) = 0.0000000001; % now take out the zeros

% normalize so in the 0-1 range
cf = cf/max(max(cf));

% plot
rdmStruct.RDM = cf;
[pearsonR, p_pearson, pointCoords] = rsa.MDSConditions_iar(rdmStruct, userOptions, localOptions);
set(gcf, 'Position',[508   462   608   484]);

import numpy as np
import rsatoolbox
import scipy.io
from scipy.io import loadmat

# This script loads in preprocessed eyetracking data for P1 and uses it to build and save out RDM files. These files are
# then reimported back into MATLAB and plotted as RDMs for figure S2c and S2d using the script figureS2cd_plotsP1.m

# Isabelle Rosenthal 2025

# note that random seed may need to be set to ensure replicable results

# directory where the output file will be saved
saveName = '../Data/stimContext_Timing_P1_RDMs_EYE.mat'

plotRDM = 0 # to plot the RDM directly in Python

# load in preprocessed data file
matFile = '../Data/stimContext_Timing_P1_EyeForRSA.mat'
matData = loadmat(matFile)

matData = matData["outStruct"]

# extract variables from the data file
segments = matData["segVec"][0][0] # was nsp x features x segments (concatenated phases/trials). Now features x segments
condLabels = matData["condLabels"][0][0][0]
condNameVec = np.array([a[0] for a in matData["condNames"][0][0][0]])
featList = np.arange(1,7)

# construct the RDM using cross-validated Mahalanobis distance (LDS) with multivariate noise normalization in every bin
rdmList = []

# Since there are different amounts of data per condition (less in realistic condition)
# iterate 1000 times, subsampling the abstract trials to match the number of realistic trials
for ind in np.arange(1000):
    # shuffle realistic trials
    Aind = np.where(condLabels==2)[0] # get all the abstract indices based on abstract ITI
    realInds = np.where(np.mod(condLabels,2)==1)[0]
    segToUse = segments[:,realInds] # take all the R data
    condLabelsToUse = condLabels[realInds]
    condNameVecToUse = condNameVec[realInds]
    np.random.shuffle(Aind)

    # subselect
    AtoUse = np.expand_dims(Aind[np.arange(sum(condLabels==1))],axis=0) # get the same number of trials as real
    AtoUse = np.concatenate((AtoUse, AtoUse+1, AtoUse+2, AtoUse+3),axis=1) # get the other phases as well as ITI
    segToUse = np.concatenate((segToUse,segments[:,AtoUse[0]]),axis=1)
    condLabelsToUse = np.concatenate((condLabelsToUse,condLabels[AtoUse[0]]))
    condNameVecToUse = np.concatenate((condNameVecToUse,condNameVec[AtoUse[0]]))

    # format the data as an rsa dataset with 6 channels
    rsaData = rsatoolbox.data.Dataset(segToUse.T,
                                      channel_descriptors= {'featList':featList},
                                      obs_descriptors={'condLabels':condLabelsToUse,'condNames': condNameVecToUse})
    # multivariate noise normalization
    noise_shrink = rsatoolbox.data.noise.prec_from_measurements(rsaData, obs_desc='condLabels', method='shrinkage_eye')

    # create rdm (representational dissimilarity matrix) with Mahalanobis distance
    rdm = rsatoolbox.rdm.calc_rdm(rsaData, method='mahalanobis', descriptor='condLabels',  noise=noise_shrink)
    outRDM = rdm.get_matrices()[0]
    rdmList.append(outRDM)

# save file out
outDict = {"RDMs": rdmList, "condLabelOrder": rdm.pattern_descriptors}
scipy.io.savemat(saveName, outDict)
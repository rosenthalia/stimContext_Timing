import numpy as np
import rsatoolbox
import scipy.io
from scipy.io import loadmat

# This script loads in preprocessed firing rate data for P2 and uses it to build and save out RDM files. These files are
# then reimported back into MATLAB and plotted as RDMs for figure 5d and 5e using the script figure5de_plotsP2.m

# Isabelle Rosenthal 2025

# note that random seed may need to be set to ensure replicable results

# directory where the output file will be saved
saveName = '../Data/stimContext_Timing_P2_RDMs.mat'

# load in preprocessed data file
matFile = '../Data/stimContext_Timing_P2_SpksForRSA.mat'
matData = loadmat(matFile)

matData = matData["outStruct"]

# extract variables from the data file
condNames = np.array([a[0] for a in matData["condNames"][0][0][0]])
chList = matData["chList"][0][0][0]
segments = matData["segVec"][0][0]

setLabels = matData["setlabelVec"][0][0]
setLabels = np.array([a[0] for a in setLabels])
condLabels = matData["phaselabelVec"][0][0]
condLabels = np.array([a[0] for a in condLabels])
condNameVec = np.array([a[0] for a in matData["phaseNamesVec"][0][0][0]])

# construct the RDM using cross-validated Mahalanobis distance (LDS) with multivariate noise normalization in every bin

# format the data as an rsa dataset with 128 channels
rsaData = rsatoolbox.data.Dataset(np.squeeze(segments[:,:].T),
                                  channel_descriptors= {'chList':np.squeeze(chList[0])},
                                  obs_descriptors={'condLabels':condLabels, 'condNames': condNameVec, 'sets':setLabels})
# multivariate noise normalization
noise_shrink = rsatoolbox.data.noise.prec_from_measurements(rsaData, obs_desc='condLabels', method='shrinkage_eye')

# create rdm (representational dissimilarity matrix) with cross-validated Mahalanobis distance
rdm = rsatoolbox.rdm.calc_rdm(rsaData, method='crossnobis', descriptor='condLabels',  #,noise=noise_shrink)
                              cv_descriptor='sets',noise=noise_shrink)
outRDM = rdm.get_matrices()[0]

# save file out
outDict = {"RDM": outRDM, "condLabelOrder": rdm.pattern_descriptors}
scipy.io.savemat(saveName, outDict)
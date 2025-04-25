README.txt

This repository contains code to reproduce analyses and figures from the paper "Visual context affects the perceived timing of tactile sensations elicited through intra-cortical microstimulation: a case study of two participants" by Rosenthal et al.

Preprocessed data to use with the code can be found at [https://doi.org/10.5281/zenodo.15284114], you will need to download the data and put it a folder called 'Data', which should be in the same directory as the 'Code' folder.

Code and data with 'P1' in the name pertains to participant P1; code and data with 'P2' in the name pertains to participant P2.

FIGURE 1: use figure1_2_S1_plots_P1.m and figure1_2_S1_plots_P2.m to plot the participant accuracies and display associated statistics.

FIGURE 2: use figure1_2_S1_plots_P1.m and figure1_2_S1_plots_P2.m to plot the rates of ICMS-elicited sensations and display associated statistics.

FIGURE 3: use figure3_plots_P1.m to plot the intensity, location, and qualitative descriptors of ICMS-elicited sensations in P1 and display associated statistics. (This data was not collected in P2.)

FIGURE 4: use figure4_plots_P1.m and figure4_plots_P2.m to plot the temporal binding window and parametric bootstrap of temporal offset data, and display associated statistics.

FIGURE 5:
 (a), (b), (c): use figure5abc_plots_P1.m and figure5abc_plots_P2.m to run the tuning analysis, plot the results, and display associated statistics.
 (d), (e): RDMs are constructed in Python using figure5de_constructRDMs_P1.py and figure5de_constructRDMs_P2.py, and the resulting save files are loaded into MATLAB for visualization using figure5de_plots_P1.m and figure5de_plots_P2.m. The analysis in Python has been prerun so if you want you can just run the MATLAB scripts directly. For python, need Python RSA toolbox (https://github.com/rsagroup/rsatoolbox). RDMS and the MDS analysis are visualized using the MATLAB rsa toolbox(https://github.com/rsagroup/rsatoolbox_matlab) which has been slightly modified for the purposes of this paper.



FIGURE S1: use figure1_2_S1_plots_P1.m and figure1_2_S1_plots_P2.m to plot the rates of ICMS-elicited sensations by timing offset and associated statistics.

FIGURE S2:
 (a), (b): use figureS2ab_plots_P1.m and figureS2ab_plots_P2.m to run the tuning analysis, plot the results, and display associated statistics for the eyetracking data.
 (c), (d): RDMs are constructed in Python using figureS2cd_constructRDMs_P1.py and figureS2cd_constructRDMs_P2.py, and the resulting save files are loaded into MATLAB for visualization using figureS2cd_plots_P1.m and figureS2cd_plots_P2.m. The analysis in Python has been prerun so if you want you can just run the MATLAB scripts directly. For python, need Python RSA toolbox (https://github.com/rsagroup/rsatoolbox). RDMS and the MDS analysis are visualized using the MATLAB rsa toolbox(https://github.com/rsagroup/rsatoolbox_matlab) which has been slightly modified for the purposes of this paper.
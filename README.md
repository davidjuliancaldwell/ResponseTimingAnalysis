### MATLAB code to produce basic figures and generate .mat files for further analysis in R and python of the stimulation response time, modified waveform, and temporal order judgement experiments in epilepsy patients. This also extracts the neural data from each type of these experiments, and contains code to perform artifact rejection and subsequently process the signals.

---

This repository contains MATLAB code to analyze the data from different sorts of stimulation response time, modified waveform, and temporal order judgement experiments in epilepsy patients

---

### Response timing and priming waveform analysis

The folder ***dataExtraction*** contains the code to extract the behavioral data from the response timing experiments for subjects 693ffd, acabb1, 2fd831, and c19968.

The main script is ***extractStimResponse_iterate.m***, which calls each individual subject script to extract and save the behavioral data.

***neuralSignalPrep*** has scripts to extract the neural data for the behavioral data once the behavioral responses have been extracted.

The main script is ***extractNeuralData_iterate.m***, which calls each individual subject script to extract the neural data for further analysis.

The folder ***neuralSignalAnalysis*** performs artifact processing, time series, and time-frequency analyses on the extracted data. This is still very much in progress.

---

### Temporal order judgement analysis

A more complete set of code for the temporal order judgement task is in the ***TOJ_shared*** project on github, but here, each subject has their individual analysis script titled ***extractNeuralData_TOJ_analysis_subjectID.m*** to extract the behavioral data subsequently saved as ***subjectID_TOJ_matlab.mat***. The subsequent neural data is extracted with ***extractNeuralData_TOJ_analyis_subjectID.m***.

***compare_TOJ_haptic.m*** is an example script comparing TOj responses to haptic only responses. 
---

### Helper functions

There are various folders with helper functions that are called by scripts in this project. These include the ***analysisFunctions***, ***visualization***, and ***behavioralDataFunctions***.

---

David J Caldwell, BSD-3 License

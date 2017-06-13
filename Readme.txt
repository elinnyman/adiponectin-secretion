---Introduction---
To be able to run the files, MATLAB and IQM tools are necessary. 
IQM tools can be aquired for free from http://www.intiquan.com/iqm-tools/

To test a simple example, run "SimplePlotExample", the plot should be equal to the one in the root folder named "Simple Plot". 
To generate the figures from the article, run "GenerateFigures". 

The data files are available in ./Plots/Data/
All model files can be found in ./Models/
All scripts (except SimplePlotExample and GenerateFigures) are available in ./Scripts/

--- Description of scripts ---
* SimplePlotExample - A script to test that everything behaves as is espected. 
* GenerateFigures - Generates the figures in the article, saves the figures in ./Plots/tmp

* CostFunction - Calculates a scalar value of the agreement between simulation and experimental data
* EstimateParameters - Tries to find the best model parameters with respect to agreement with data, for set of model structures
* Intialization - Sets variables necessary for other scripts
* MergeValidParameters - Merges multiple tables of parameter sets
* PlotFigures - Plots multiple figures, such as agreement with data, behavior of the input etc
* PlotIntervals - Plots the profile-likelihood figures. 
* PlotUncertainity - Plots the the agreement with data, for multiple parametersets
* QuantifyPrediction - Quantifies the uncertainity in the models prediction of new data. 
* RecheckParameters - Reruns a collection of parameter sets, to dubblecheck that all parameter sets have a sufficient agreement with data.
* SimulateExperiments - Simulates a set of experiment with a model, to recreate the experiments performed when collectning the experimental data

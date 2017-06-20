% The scripts makes sure that the scripts from the right folders, and sets
% necessary parameters. Needs to have "modelName" defined before it
% can be run. 

%% Makes sure the scripts are run from the right folder.
if exist('../../Scripts/SimulateExperiments.m','file')
    cd ../..
elseif exist('../Scripts/SimulateExperiments.m','file')
    cd ..
elseif exist('../../../Scripts/SimulateExperiments.m','file')
    cd ../../..
elseif ~exist('./Scripts/SimulateExperiments.m','file')
    error('Running from wrong folder, please change to root folder.')
    return
end

if ~exist('IQMmakeMEXmodel.m','file')
    error('IQM tools missing, please install from http://www.intiquan.com/iqm-tools/.')
    return
end
%% Sets the necessary parameters to globals
global model
global expData
global highCaExpData
global time
global highCaTime
global expData_cAMP
global time_cAMP
global validParameters
global dgf
global lb
global ub
global wTotalData

%%
clear mex
addpath('Data')
addpath(genpath([pwd '/Models']))
addpath('Scripts')

load expData
load highCaExpData
load expData_cAMP

if ~exist(['./Models/' modelName],'dir')
    mkdir(['./Models/' modelName])
end

cd Models
IQMmakeMEXmodel(IQMmodel([modelName '.txt'])) % Compiles a MEX file of the model. 
cd(modelName)

if ~exist('./Plots','dir')
    mkdir('./Plots')
end
if ~exist('./Plots/tmp','dir')
    mkdir('./Plots/tmp')
end
if ~exist('./Results','dir')
    mkdir('./Results')
end

model=str2func(modelName); % Sets the model as a function.
pNames=IQMparameters(model); % Gets the parameter names from the model. 

if ~exist('wTotalData','var') || isempty(wTotalData) % Sets the choice of datasets to the total, if no other choice is made prior
 wTotalData=1;
end

if ~exist('d','var') || isempty(d) % Sets the choice of datasets to the total, if no other choice is made prior
 d=0.1;
end

if ~exist('optParam','var') % Uses the previous best solution as start guess. If no solution is available, uses the parmetervalues from the model file
    [pNames, startGuess]=IQMparameters(model);
    startGuess=startGuess(1:end-4)';
else
    startGuess=optParam;
end

dgf=numel(expData.meanValues)+ (numel(highCaTime)+numel(time_cAMP))*wTotalData-1; % Sets the degrees of freedom for the chi2-limit. Number of measured time points in the datset used -1 scaling parameter. 
nParams=length(startGuess);
ub=log(1e4)*ones(1,nParams-6); % Sets the upper bounds of the parameter values (in log-space). 
lb=-ub; % Sets the lower bounds.

ub=[ub log([0.7    1+d 1+d 0.23  1     0.1 ])]; % Add upperbounds for the rate corresponding to the half-maximal time, the difference in time between different stimuli, and basal concentrations of Ca, ATP, cAMP
lb=[lb log([0.0035 1-d 1-d 0.035 1e-3  1e-6])]; % Add lower bounds 

% uncomment to limit cAMP to experimentally determined values
% lb=[lb log([0.0035 1-d 1-d 0.035 1e-3  0.00198823])]; % Add lower bounds, with limitation on cAMP

chi=chi2inv(0.95,dgf); % Sets the chi2-limit

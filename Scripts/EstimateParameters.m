% This script estimates the parameters for a set of models.

clear all
global validParameters


models={ %List of models to have parameters estimated.
'ATP_Sat_Endo'
'OriginalModel';
'ATPdep';
% 'Sat';
% 'Endo';
% 'ATP_Sat';
% 'ATP_Endo';
% 'Sat_Endo';
};
nRestarts=10; % Determines how many times the optimization algorithm will be run. 
weights=[0; 1]; % Determines which datasets will be used, 0 = original, 1 = total.
d=0.1; % Sets how much difference is allowed in the time to half maximal response. 0.1 corresponds to 10% deviation. 

cd Models
for i =1:length(models) % Compiles all the models into MEX files.
    IQMmakeMEXmodel(IQMmodel([models{i} '.txt'])) 
end
cd ..
 
psOpts = optimoptions(@particleswarm,'Display','iter'); % Particle swarm optimization settings
saOpts = optimoptions(@simulannealbnd,'HybridFcn',@fmincon,'Display','iter'); % Simulated annealing optimization settings

nWeights=size(weights,1);
results=table(reshape(repmat(models',nWeights,1),1,length(models)*nWeights)',repmat(weights,length(models),1), zeros(length(models)*size(weights,1),nRestarts),'VariableNames',{'Model','Weights','Costs'}); % Contructs a table, containing a combination of models and weights, and costs
for i=1:height(results)
    clear mex
    clear optParam
    validParameters=[]; % Clears the table of parameter sets passing a chi2 test
    modelName=results.Model{i}; %Selects with models to
    
    wTotalData=results.Weights(i,1);
    Initialization % Sets necessary variables, depending on which model is being run. 
    
    for j=1:nRestarts
        [optParamPS, minfunPS]=particleswarm(@CostFunction, nParams, lb,ub,psOpts); % Optimize parameters using particle swarm.
        [optParam, minfun]=simulannealbnd(@CostFunction, optParamPS, lb,ub,saOpts); % Using the best result from particle swarm as an intial guess, and runs an simulated annealing optimization algorithm. 
        results.Costs(i,j)=minfun; 
        save(sprintf('./Results/opt(%.2f), %s.mat',minfun,  datestr(now,'yymmdd-HHMMSS')),'optParam') % Saves the best parameters found.
        if ~isempty(validParameters)
            save(sprintf('./Results/valid, %s.mat',  datestr(now,'yymmdd-HHMMSS')),'validParameters') % Saves all parameters passing a chi2 test, if any exist. 
        end
        
    end
end

save(sprintf('../Results %s.mat', datestr(now,'yymmdd-HHMMSS')),'results') % Saves the results table
writetable(results, sprintf('../Results %s.xlsx', datestr(now,'yymmdd-HHMMSS'))) % Save the results table as a .xslx file. 
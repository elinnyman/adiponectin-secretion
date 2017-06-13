function [ cost, endoExo ] = CostFunction( param )
%This function is used to give a measure of how well the simulation agrees
%with the experimental data. Estimated with the sum of squared residuals,
%normalized with the SEM. 

global expData
global highCaExpData
global expData_cAMP
global time
global highCaTime
global time_cAMP
global validParameters
global dgf
global wTotalData


originalCost=0;
cAMPCost=0;

try
    [sim,~,endoExo, scale]=SimulateExperiments(param, 0:1:720); % Simulates the experiments.
    
    %% Compares the original experimental data with the corresponding simulations
    tInd=ismember(sim.measures(end,:),time); % Finds which time points to use when comparing simulation to experimental data
    for i=1:height(expData)
        originalCost=originalCost+nansum((sim.measures(i,tInd)-expData.meanValues(i,:)).^2./expData.SEMValues(i,:).^2); % sum( (simulation - experimental)^2/SEM^2)
    end
    
    %% Compares the new high Ca experimental data with the corresponding simulation

    i=4+1;
    tInd=ismember(sim.measures(end,:),highCaTime);
    CaCost=nansum((sim.measures(i,tInd)-highCaExpData.meanValues).^2./highCaExpData.SEMValues.^2); % calculates the cost for 
    
    %% Compares the 0 cAMP experimental data with the corresponding simulation

    tInd=ismember(sim.measures(end,:),time_cAMP);
    cAMPCost=cAMPCost+nansum((sim.measures(1+i,tInd)-expData_cAMP.meanValues(1,:)).^2./expData_cAMP.SEMValues(1,:).^2);
    cost=originalCost+(CaCost+cAMPCost)*wTotalData;
    
    %% Saves parameters passing a chi2 test (shifted by +15 to get a sense of parameters close to the boundry, later removed)
    if cost<chi2inv(0.95,dgf)+15
        validParameters=[validParameters; table(param, cost, originalCost, CaCost, cAMPCost, endoExo, scale, 'VariableNames', {'Parameters', 'Cost','OriginalCost', 'CaCost','cAMPCost', 'ExoEndo','scale'})];
    end
catch err
    disp(getReport(err))
    cost=inf;
end
end


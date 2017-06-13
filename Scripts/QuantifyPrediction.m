% This script quantifies the uncertainity in the simulations of released
% adiponectin. 

vp=validParameters;
vp(vp.Cost>chi2inv(0.95,dgf-size(vp.Parameters,2)),:)=[]; % removes all parameter set where the cost is larger than the chi2 limit.

change=[1            1 1; %Defines how much the inputs should increase given different stimuli (as fold change over basal).
        1.5*1.25     1 6.5*1.25;
        1.5*0.75     1 6.5*0.75;
        0            1 6.5*1.25;
        0            1 6.5*0.75;
        10           1 6.5*1.25;
        1.5+1.5*0.75 1 6.5*0.75;
    ];

design={'Control';'fsk+/+'; 'fsk-/-' ; 'noCa+'; 'noCa-';'Iono10+'; 'Iono1.5-';'time'};

values=nan(height(vp),size(change,1)-1);
for i=1:height(vp) % Simulates the adiponectin release for all parameter sets.
    param=vp.Parameters(i,:);
    basalValues=exp(param(end-2:end));
    parameters=repmat(basalValues,size(change,1),1).*change;
    [simulatedExperiments, allSimData]=SimulateExperiments(param, [0 time 30*60], design, parameters);
    simVals=[];
    for k=1:length(allSimData)
        simVals=[simVals allSimData(k).statevalues(end,ismember(allSimData(k).states,'Adiponectin'))];
    end
    simVals=simVals./simVals(1); % Normalize with time 0. 
    values(i,:)=simVals(2:end);
    
    if mod(i,1000)==0
        fprintf('Done with %i of %i \n', i, height(vp))
    end
end

save values values % saves the values simulated as backup. 
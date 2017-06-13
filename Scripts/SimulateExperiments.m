function [ simulatedExperiments, allSimData, endoExo, scale] = SimulateExperiments(param, times, design,parameters)
% This function simulates the experiments done. It takes a parameter set,
% requested time points, names of the experimental design and experimental
% parameters as inputs. 
% Returns the simulated timeseries, a table of all simulation structures,
% endo-/exocytosis and the scaling parameter. 
%% Initial setups
global model
global expData
global time


if nargin==1
    times=time;
end

if times(1)~=0
    times=[0 times];
end

endoExo=zeros(1,6);
param=exp(param);


if nargin==4
    parameters=[repmat(param,size(parameters,1),1) parameters];
elseif nargin==3
    disp('Design must be paired with parameters! Using default.')
elseif nargin <3
    design={'Ca+ATP';'Ca+noATP';'noCa+ATP';'noCa+noATP';'10xCa+ATP';...
        '1.5 Ca+ATP+nocAMP';'noCa+ATP+nocAMP';'noCa+noATP+nocAMP';'time'};
    parameters=[1.5   3  0.1;
        1.5   0  0.1;
        0  3  0.1;
        0  0  0.1;
        15    3  0.1;
        1.5   3  0   ; 
        0  3  0   ;
        0  0  0   ];
    
    parameters=[repmat(param,size(parameters,1),1) parameters];
    
end
parameters(~cellfun(@isempty,strfind(design,'noCa')),[end-5, end-2])=param(end-2);

measures=[];
simulations=[];
timeInd=ismember(0:1:times(end),times);

%% Simulate experiments
for i=1:size(parameters,1)
    sim=model(0:1:times(end),[],[parameters(i,:) 1]); %Simulate the model
    maxT=find(sim.variablevalues(:,1)==max(sim.variablevalues(:,1)),1); %Find the time and collect values for endo-/exocytosis contribution
    if strcmp(design{i},'Ca+ATP')
        endoExo(1:2)=sim.variablevalues(maxT,2:3);
    elseif strcmp(design{i},'Ca+noATP')
        endoExo(3:4)=sim.variablevalues(maxT,2:3);
    elseif strcmp(design{i},'10xCa+ATP')
        endoExo(5:6)=sim.variablevalues(maxT,2:3);
    end
    
    sim.time=sim.time(timeInd);
    sim.variablevalues=sim.variablevalues(timeInd,:);
    sim.statevalues=sim.statevalues(timeInd,:);
    simulations=[simulations sim]; % Collect the full simulation structure, but only for the time points requested. 
    
    if length(times)==length(time)+1 % Stores the simulated values in a matrix. 
        measures=[measures simulations(end).variablevalues(2:end,1)];
        simulations(end).time(1)=[];
    else
        measures=[measures simulations(end).variablevalues(:,1)];
    end
    
    
end

if size(simulations,2)>3 % scales the simulation against the experimental data. 
    meanValues=expData.meanValues';
    sims=measures(ismember(simulations(end).time,time),1:height(expData));
    
    sims=reshape(sims,1,numel(sims))';
    meanValues=reshape(meanValues,1,numel(meanValues))';
    
    sims(isnan(meanValues))=[];
    meanValues(isnan(meanValues))=[];
    
    scale=lscov(sims,meanValues);
    if scale<0
        scale=1;
    end
else
    scale=1;
end
measures=[measures'.*scale; simulations(end).time];

simulatedExperiments=table(design, measures);
allSimData=simulations;
end




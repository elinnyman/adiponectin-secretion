% Recheck validParameters
% This script loads a set of saved parameters, then simulates and checks
% again if the parameterset can passa chi2-test. All parametersets that
% pass the test will be saved in a new .mat file. 


global validParameters;
modelName='ATP_Sat_Endo';
d=0.1;
Initialization

%% Load and remove invalid parameter sets. 
load('C:\Users\willo18\Documents\Projects\Adiponectin\Models\ATP_Sat_Endo\Results\validParameters, all6.mat')
valid=validParameters;
validParameters1=unique(valid,'rows');

invalid=find(sum(validParameters1.Parameters<repmat(lb,size(validParameters1,1),1),2)>0);
invalid=[invalid; find(sum(validParameters1.Parameters>repmat(ub,size(validParameters1,1),1),2)>0)];
invalid=unique(invalid);
validParameters1(invalid,:)=[];
validParameters=[];

%% Rerun all valid parameter sets
for i=1:height(validParameters1)
    cost=CostFunction(validParameters1.Parameters(i,:));
    if mod(i,50000)==0 && ~isempty(validParameters) % Saves the table of parameter sets every 50 000 iteration. 
        disp(['done with ' num2str(i) ' of ' num2str(height(validParameters1)) ' - ' num2str(cost)]);
        save(sprintf('./Results/valid, %s.mat',  datestr(now,'yymmdd-HHMMSS')),'validParameters')
        validParameters=[];
    end
end

if  ~isempty(validParameters)
    save(sprintf('./Results/valid, %s.mat',  datestr(now,'yymmdd-HHMMSS')),'validParameters')
end

function[valid]= MergeValidParameters(path)
% This function is used to merge all different "validParameters" files,
% into one big matrix. Takes a path to the folder containing the files as
% input.

addpath(path)
files=dir([path '\*.mat']); % Finds files.
valid=[];
for file=files'
    load(file.name);
    valid=[valid; validParameters;]; % Stacks all tables into one large table.
    valid=sortrows(valid,2,'ascend'); % Sorts based on the agreement with data (Cost)´
end
valid=unique(valid,'rows'); % Removes duplicate rows.
valid=sortrows(valid,2,'ascend'); % Sorts based on the agreement with data again. 
end
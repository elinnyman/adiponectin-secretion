if exist('IQMmakeMexModel','file') 
    modelName='ATP_Sat_Endo';
    load('./Models/ATP_Sat_Endo/Results/opt (24.91).mat')
    cd Models
    IQMmakeMEXmodel(IQMmodel([modelName '.txt']))
    cd ..
    
    d=0.1;
    wHi=1;
    wcAMP=1;
    
    addpath([pwd '/Data'])
    addpath(genpath([pwd '/Models']))
    addpath([pwd '/Scripts'])
    
    Initialization
    
    PlotFigures(optParam)
    cd ../..
else
    disp('IQMtoolbox might be necessary to run all scripts. Please install it first.')
    disp('It can be freely downloaded from: http://www.intiquan.com/iqm-tools/')
end
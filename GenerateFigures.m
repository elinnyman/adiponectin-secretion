warning off
clear all
clear mex

%% Initial Setup. 

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
home=pwd;

if ~exist('IQMmakeMEXmodel.m','file')
    error('IQM tools missing, please install from http://www.intiquan.com/iqm-tools/.')
    return
end
addpath('Data')
addpath(genpath([pwd '/Models']))
addpath('Scripts')
d=0.1;

%% Figure 2 Original+ATPdep
wTotalData=0;
modelName='OriginalModel';
load('./Models/OriginalModel/Results/opt(62.13), 170503-202614.mat')
Initialization;
cd(home)
PlotFigures( optParam, 1, 1, 1, 0, 0, 0 )

modelName='ATPdep';
load('./Models/ATPdep/Results/opt(30.57), 170503-205819.mat')
Initialization;
cd(home)
PlotFigures( optParam, 1, 1, 1, 0, 0, 1 )

%% Setup for final model figures 
close all
modelName='ATP_Sat_Endo';
load('./Models/ATP_Sat_Endo/Results/opt (10.59).mat')
load('./Models/ATP_Sat_Endo/Results/validParameters, all.mat');
wTotalData=1;
Initialization;
cd(home)

%% Figure 3 ATP_Sat_Endo
PlotUncertainity(validParameters)


%% Figure 4 - patch-clamp vs extracellular stimuli
figure(4)
meanValue= 2.050012250e-14 * 1e15;
SEMValue=3.509235835e-15 * 1e15;

bar(2,expData.meanValues(1,2),'FaceColor',[0.6 0.6 0.6])
hold on
errorbar(2,expData.meanValues(1,2),expData.SEMValues(1,2),'ko','MarkerFace','auto','linewidth',2)
bar(1,meanValue,'FaceColor',[0.3 0.3 0.3])
errorbar(1,meanValue,SEMValue,'ko','MarkerFace','auto','linewidth',2)
set(gca, 'XTick', 1:2, 'XTickLabel', {'FSK/IBMX','1.5 µM Ca^{2+}'});
ylabel('\DeltaC/\Deltat (fF/s)')
box off
set(gca,'FontName', 'Arial')
set(gca,'FontSize',15)
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0.1,0.1,0.2,0.4]);  


savefig('./Plots/tmp/Figure 4, bars')
print('./Plots./tmp/Figure 4, bars.pdf','-dpdf')
print('./Plots/tmp/Figure 4, bars.png','-dpng')


%% Figure 5 Prediction/validation
% QuantifyPrediction; %Used to collected the prediction values in
% "values.mat"
load('./Models/ATP_Sat_Endo/values.mat')

vp=validParameters;
vp.Cost=validParameters.OriginalCost+validParameters.CaCost+validParameters.cAMPCost;
vp(vp.Cost>chi2inv(0.95,dgf-size(vp.Parameters,2)),:)=[];
ind = any(vp.Parameters>repmat(ub,height(vp),1),2) | any(vp.Parameters<repmat(lb,height(vp),1),2) | exp(vp.Parameters(:,end))>0.1 | exp(vp.Parameters(:,end))<0.00198823;
vp(ind,:)=[];
values(ind,:)=[];
change=[1 1 1;
    1.5*1.25 1 6.5*1.25;
    1.5*0.75 1 6.5*0.75;
    0 1 6.5*1.25;
    0 1 6.5*0.75;
    10 1 6.5*1.25;
    1.5+1.5*0.75 1 6.5*0.75;
    ];
design={'Control';'fsk+/+'; 'fsk-/-' ; 'noCa+'; 'noCa-';'Iono10+'; 'Iono1.5-';'time'};

PlotPrediction;

%% Figure 6, Endocytosis contribution intervals
close all
vp=validParameters;
vp.Cost=validParameters.OriginalCost+validParameters.CaCost+validParameters.cAMPCost;
vp(vp.Cost>chi2inv(0.95,dgf-size(vp.Parameters,2))+10,:)=[];
ind = any(vp.Parameters>repmat(ub,height(vp),1),2) | any(vp.Parameters<repmat(lb,height(vp),1),2) | exp(vp.Parameters(:,end))>0.1 | exp(vp.Parameters(:,end))<0.00198823;
vp(ind,:)=[];
PlotInterval(vp,4)


%% Figure 6, quantification example
endo=vp.ExoEndo(:,1:2:5);
exo=vp.ExoEndo(:,2:2:6);
vp.ExoEndo=endo./(endo+exo);
vp=sortrows(vp,'Cost');

ind=find(vp.ExoEndo(:,1)<0.10,1);

param=vp.Parameters(ind(1),:);
PlotFigures(param, 0, 0, 0, 0, 1, 0)

warning on
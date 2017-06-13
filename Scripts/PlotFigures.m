function [ ] = PlotFigures( optParam, data, states, input, steady, measure, overlay )
% This function plots multiple different and relevant figures.

global time
global model
global scale

if nargin==1 % Sets standard values of which plots to generate. Note: Will fail if giving more than one, less than all inputs.
    data=true;
    states=false;
    input=false;
    steady=false;
    measure=false;
    overlay=false;
end

if ~overlay % Closes all plots if wanted. Otherwise keeps the elements in the previous plots (useful if plotting and comparing two models)
    close all
end
CostFunction(optParam)
design={'1.5 Ca+ATP';'1.5 Ca+noATP';'noCa+ATP';'noCa+noATP';'15 Ca+ATP';'0.7 Ca+ATP';'0.25 Ca+ATP';'150 Ca+ATP';'10 Ca+ATP';'1.5 Ca+ATP+nocAMP';'noCa+ATP+nocAMP';'noCa+noATP+nocAMP';'time'};

parameters=[1.5   3  0.1;
    1.5   0  0.1;
    0  3  0.1;
    0  0  0.1;
    15    3  0.1;
    0.7    3  0.1;
    0.25    3  0.1;
    150    3  0.1;
    10    3  0.1;
    1.5   3  0   ;
    0  3  0   ;
    0  0  0   ];

[simulatedExperiments, allSimData, ~, scale]=SimulateExperiments(optParam, 0:1:time(end), design, parameters);

if data
    PlotData(simulatedExperiments, states)
end
if states
    PlotStates(allSimData,simulatedExperiments.design)
end
if input
    if strfind(func2str(model),'Can')
        PlotInput(allSimData, simulatedExperiments.design, optParam(end))
    else
        PlotInput(allSimData([1:4 9:12]), simulatedExperiments.design([1:4 9:12 ]))
    end
end

if steady
    PlotSteadyState(optParam)
end
if measure
    PlotEndocytosis(allSimData)
end

end

function [] = PlotData(sim, states)
%Produce a simple plot of the agreement with data, similiar to the one in the article.
global model

figure(1);
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);
hold on
global time
global expData
global highCaExpData
global highCaTime
global expData_cAMP
global time_cAMP

if nargin == 0
    sim=[];
end
simtime=sim.measures(end,:)./60;

if states
    n=3;
else
    n=2;
end

for i=1:height(expData)
    subplot(n,2,i)
    hold on
    if strcmp(func2str(model),'ATPdep')
        plot(simtime,sim.measures(i,:),'--','color','b', 'LineWidth', 3)
    else
        plot(simtime,sim.measures(i,:),'color',[0.7,0.7,0.7], 'LineWidth', 3)
    end
    hold on
    errorbar(time./60,expData.meanValues(i,:),expData.SEMValues(i,:),'ko','MarkerFaceColor','k')
    title(sim.design{i})
    axis([0 12.5 -1 24])
    set(gca,'ytick',[0 10 20], 'xtick',0:2:12)
    xlabel('Time (Min)')
    
    ylabel('\DeltaC/\Deltat (fF/s)')
    box off
    set(gca,'FontName', 'Arial')
    %
end

%% Plot high Ca Data
figure(6)
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);
hold on
if  strcmp(func2str(model),'ATPdep')
    plot(simtime,sim.measures(5,:),'b--', 'LineWidth', 3)
elseif strcmp(func2str(model),'ATP_Ca')
    plot(simtime,sim.measures(5,:),':','color',[0.7,0.7,0.7], 'LineWidth', 3)
else
    plot(simtime,sim.measures(5,:),'color',[0.7,0.7,0.7], 'LineWidth', 3)
end
errorbar(highCaTime./60,highCaExpData.meanValues,highCaExpData.SEMValues,'ko','MarkerFaceColor','k')
title(sim.design{5})
xlabel('Time (Min)')
ylabel('\DeltaC/\Deltat (fF/s)')
box off
set(gca,'FontName', 'Arial')
axis([0 8.5 -1 24])


%%
figure(28)
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);

plot(sim.measures(end,:)./60, sim.measures(10,:)./60,'linewidth',2,'color','k')
hold on
plot(sim.measures(end,:)./60, sim.measures(11,:)./60,'linewidth',2,'color','r')
plot(sim.measures(end,:)./60, sim.measures(12,:)./60,'b--','linewidth',2)
errorbar((time_cAMP-5)./60,expData_cAMP.meanValues(1,:),expData_cAMP.SEMValues(1,:),'ks','markerfacecolor','auto','linewidth',2)
errorbar(time_cAMP./60,expData_cAMP.meanValues(2,:),expData_cAMP.SEMValues(2,:),'ro','markerfacecolor','auto','linewidth',2)
errorbar((time_cAMP+5)./60,expData_cAMP.meanValues(3,:),expData_cAMP.SEMValues(3,:),'b^','markerfacecolor','auto','linewidth',2)

xlabel('Time (Min)')
ylabel('\DeltaC/\Deltat (fF/s)')
axis([0 12.5 -2.1 2.1])
box off
legend(sim.design(10:12));

%% save plots
figure(1)
set(findall(gcf,'-property','FontSize'),'FontSize',15)
savefig ./Plots/tmp/data.fig
print -dpdf -r0 ./Plots/tmp/data.pdf
print -djpeg -r0 ./Plots/tmp/data.jpg

figure(6)
set(findall(gcf,'-property','FontSize'),'FontSize',15)
savefig ./Plots/tmp/data_highCa.fig
print -dpdf -r0 ./Plots/tmp/data_highCa.pdf
print -djpeg -r0 ./Plots/tmp/data_highCa.jpg
end

function [] = PlotStates(sims, design)
% Plots the releable pool (Fig 2E,F)
global model

figure(1);
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);
hold on


state=5;
for experiment=1:2
    subplot(3,2,4+experiment)
    if strcmp(func2str(model),'ATPdep')
        plot(sims(experiment).time./60,sims(experiment).statevalues(:,state)./sims(experiment).statevalues(1,state)*100,'--','color','b','linewidth',3)
    elseif strcmp(func2str(model),'ATP_Ca')
        plot(sims(experiment).time./60,sims(experiment).statevalues(:,state)./sims(experiment).statevalues(1,state)*100,':','color',[0.7,0.7,0.7],'linewidth',3)
    else
        plot(sims(experiment).time./60,sims(experiment).statevalues(:,state)./sims(experiment).statevalues(1,state)*100,'color',[0.7,0.7,0.7],'linewidth',3)
    end
    box off
    axis([0 12 0 110])
    hold on
    
    title([sims(experiment).states{state} ' - ' design{experiment}])
    
    ylabel('Releaseable pool, % of Max')
    xlabel('Time (Min)')
    set(gca,'ytick',0:50:100, 'xtick',0:2:12)
end
figure(1)
savefig ./Plots/tmp/data.fig
print -dpdf -r0 ./Plots/tmp/data.pdf
print -djpeg -r0 ./Plots/tmp/data.jpg
end

function [] = PlotInput(sims, design)
% Plots the time behavior of the inputs

figure;

set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);

if nargin<2
    design={'Ca+ATP';'Ca+noATP';'noCa+ATP';'noCA+noATP'};
end

pos=1;
inputs=3;
nExperiments=length(design);

for experiment=1:nExperiments
    for state=1:inputs
        subplot(nExperiments,inputs,pos)
        plot(sims(experiment).time./60,sims(experiment).statevalues(:,state),'linewidth',3)
        box off
        axis tight
        
        if experiment==1
            title(sims(experiment).states{state})
        end
        if state==5
            ylabel(design{experiment})
        end
        
        set(gca,'fontsize',12)
        pos=pos+1;
    end
    
    set(gca,'fontsize',12)
end
figure(3)
savefig ./Plots/tmp/input.fig
print -dpdf -r0 ./Plots/tmp/input.pdf
print -djpeg -r0 ./Plots/tmp/input.jpg
end

function[]=PlotSteadyState(optParam)
% Plots a simulation run for a long time, to check steady state. 
global model

figure(9)
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);
sim=model(0:1:400,[],[exp(optParam) 1.5 3 0.1 0]);
n=length(sim.states);
a=ceil(sqrt(n));
b=ceil(n/a);
for j=1:n
    subplot(a,b,j)
    plot(sim.time,sim.statevalues(:,j));
    title(sim.states(j));
end
savefig ./Plots/tmp/steady.fig
print -dpdf -r0 ./Plots/tmp/steady.pdf
print -djpeg -r0 ./Plots/tmp/steady.jpg
end

function [] = PlotEndocytosis(sims)
% Plots the Recorded signal, and the contributions from the endocytosis and
% the exocytosis.
figure(5);
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);

global scale
hold on
design={'Ca+ATP';'Ca+noATP';'noCa+ATP';'noCA+noATP';'highCa+ATP';'time'};

experiments=[1 2 5];
vars=[1 2 3];
for i=1:length(experiments)
    sims(experiments(i)).variablevalues(:,2)=-sims(experiments(i)).variablevalues(:,2);
    subplot(3,2,experiments(i))
    for measures=vars
        plot(sims(experiments(i)).time./60,sims(experiments(i)).variablevalues(:,measures).*scale,'linewidth',3)
        
        hold on
        
    end
    title(design{experiments(i)})
    axis tight
    set(gca,'ytick',[0 10 20], 'xtick',0:2:12)
    xlabel('Time (Min)')
    
    ylabel('\DeltaC/\Deltat (fF/s)')
    box off
    set(gca,'FontName', 'Arial')
    legend({'Recorded signal','Endocytosis','Exocytosis'})
end


figure(5)
set(findall(gcf,'-property','FontSize'),'FontSize',15)

savefig ./Plots/tmp/endocytosis.fig
print -dpdf -r0 ./Plots/tmp/endocytosis.pdf
print -djpeg -r0 ./Plots/tmp/endocytosis.jpg
end


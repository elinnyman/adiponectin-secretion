
%% Plot prediction
% Plots the prediction, with the values from QuanitifyPrediction"
close all

[~, ind1]=min(values);
[~, ind2]=max(values);
ind=reshape([ind2(1:2:5); ind1(2:2:6)],1,length(ind1));

minVal=min(values(:,2:2:6)); % 'fsk+'; 'fsk-' ; 'noCa+'; 'noCa-';'Iono10+'; 'Iono1.5-'
maxVal=max(values(:,1:2:5));

param=vp.Parameters(ind,:);

%% Plot B-C
figure;
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);

subplot(3,2,3)
h=fill([0.75 1.25 1.25 0.75]',[0.975 0.975 1.025 1.025]','k');
set(h,'facealpha',1,'EdgeColor','None')
hold on
h=fill([1.75 2.25 2.25 1.75]',[minVal(1) minVal(1) maxVal(1) maxVal(1)]',[255 120 0]./256);
set(h,'linewidth',3)
h=fill([2.75 3.25 3.25 2.75]',[minVal(2) minVal(2) maxVal(2) maxVal(2)]',[255 120 0]./346);
set(h,'linewidth',3)
h=fill([3.75 4.25 4.25 3.75]',[minVal(3) minVal(3) maxVal(3) maxVal(3)]',[255 120 0]./436);
set(h,'linewidth',3)

title('Model Prediction')
ylabel('Released Adiponectin\newline(fold change)')
axis([0.5 4.5 0 10])
box off
set(gca, 'XTick', 1:4, 'XTickLabel', {'Control', 'FSK/IBMX', 'FSK/IBMX+BAPTA','FSK/IBMX+Iono'},'XTickLabelRotation',45);
set(gca, 'Position', get(gca, 'OuterPosition') - ...
     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

subplot(3,2,4)
errorbar(1:4, [1 3.2 3.21 5.9], [0 0.6 0.5 0.9],'ko','markerface','auto','linewidth', 2)
hold on
title('Experimental Data')
ylabel('Released Adiponectin\newline(fold change)')
set(gca, 'XTick', 1:4, 'XTickLabel', {'Control', 'FSK/IBMX', 'FSK/IBMX+BAPTA','FSK/IBMX+Iono'},'XTickLabelRotation',45);
axis([0.5 4.5 0 10])
set(gca, 'Position', get(gca, 'OuterPosition') - ...
     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%% Plot D-E
simTime=unique([(0:1:60)*60 time]);
t30=ismember(simTime, 30*60);
stateInd=ismember(IQMstates(model),'Adiponectin');
n=length(design)-1;
adiTS=table(nan(2,length(simTime)),nan(2,length(simTime)),'VariableNames',{'Control','FSK'});


for i=1:size(param,1) % Simulates the adiponectin release for all parameter sets.
    basalValues=exp(param(i,end-2:end));
    parameters=repmat(basalValues,size(change,1),1).*change;
    [~, allSimData]=SimulateExperiments(param(i,:), simTime, design([1:3 end]), [parameters(1:3,:)]);
    
    adi=[allSimData(1).statevalues(:,stateInd)'; allSimData(2).statevalues(:,stateInd)'; allSimData(3).statevalues(:,stateInd)']./allSimData(1).statevalues(t30,stateInd);
    adiTS{1,1}=max([adiTS{1,1}; adi(1,:)]);
    adiTS{2,1}=min([adiTS{2,1}; adi(1,:)]);
    adiTS{1,2}=max([adiTS{1,2}; adi(2:3,:)]);
    adiTS{2,2}=min([adiTS{2,2}; adi(2:3,:)]);
end


data=table([0 1.41301E-05	9.69864E-06	2.32011E-05	2.23375E-05; 0 2.42616E-05	2.58242E-05	3.91642E-05	5.65937E-05]*1e6, [0 2.18628E-06	6.68114E-07	1.83264E-06	1.8591E-06; 0 6.32098E-06	1.51956E-06	3.01738E-06	4.99207E-06]*1e6,'variablenames',{'MeanValues','SEMValues'});
data{:,:}=data{:,:}/data.MeanValues(1,3);
data.SEMValues(1,3)=0;

alpha=0.25;


hold on
box off
col=[1 1 1]; % black
x=simTime/60;
subplot(3,2,5)
y1=adiTS.Control(1,:);                      %create first curve
y2=adiTS.Control(2,:);                  %create second curve

X=[x,fliplr(x)];                %create continuous x value array for plotting
Y=[y1,fliplr(y2)];              %create y values for out and then back
h1=fill(X,Y,col);                  %plot filled area
 set(h1, 'EdgeColor','None')
hold on
plot(x,y1,'k','linewidth',1.5)
plot(x,y2,'k','linewidth',1.5)
title('Model Prediction')
ylabel('Released Adiponectin\newline(fold over control at 30 min)')
xlabel('Time (min)')
set(gca, 'XTick', 0:15:60);
box off
axis([0 65 0 10])
 set(gca, 'Position', get(gca, 'OuterPosition') - ...
      get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

col=[255 120 0]./256; % orange
y1=adiTS.FSK(1,:);       %create first curve
y2=adiTS.FSK(2,:);
X=[x,fliplr(x)];                %create continuous x value array for plotting
Y=[y1,fliplr(y2)];              %create y values for out and then back
h1=fill(X,Y,col);                  %plot filled area
set(h1, 'EdgeColor','None')
hold on

x=[0 15 30 45 60];                  %initialize x array
    subplot(3,2,6)
errorbar(x, data.MeanValues(1,:), data.SEMValues(1,:),'ko','markerfacecolor','auto','linewidth',1.5)
hold on
errorbar(x, data.MeanValues(2,:), data.SEMValues(2,:),'k^','color',col,'markerface','auto','linewidth',1.5)
set(gca,'FontName', 'Arial')
title('Experimental Data')
ylabel('Released Adiponectin\newline(fold over control at 30 min)')
xlabel('Time (min)')
set(gca, 'XTick', 0:15:60);
box off
axis([0 65 0 10])
 set(gca, 'Position', get(gca, 'OuterPosition') - ...
      get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

%% Plot A
colors=reshape([copper(3), copper(3)]',3,6)';
adipoInd=ismember(IQMstates(model),'Adiponectin');
simValues=[];

for k=1:size(param,1)
   basalValues=exp(param(:,end-2:end));
parameters=basalValues.*change(2:end,:);

    [simulatedExperiments, allSimData]=SimulateExperiments(param(k,:), 0:1:30*60, design([1; k+1;end]), [basalValues(k,:); parameters(k,:)]);   
    
    simValues=[simValues; allSimData(2).statevalues(:,adipoInd)'./allSimData(1).statevalues(end,adipoInd)'];
    hold on
end

titles={'FSK/IBMX'; 'BAPTA'; 'Ionomycin'};
for k=1:3
    subplot(3,3,k)
    h=fill([allSimData(1).time fliplr(allSimData(1).time)]/60, [simValues(2*k-1,:) fliplr(simValues(2*k,:))],[255 120 0]./(256+90*(k-1)));
    set(h,'facealpha',1,'EdgeColor','None')
    hold on
    errorbar(30,mean(simValues(2*k-1:2*k,end)), (simValues(2*k-1,end)-simValues(2*k,end))/2,'k','markerface','auto','linewidth',2);
    title(titles{k})
    ylabel('Released Adiponectin\newline(fold change)')
    xlabel('Time (min)')
    axis([0 32 0 10])
    
    box off
        set(gca,'FontName', 'Arial')
   set(gca, 'Position', get(gca, 'OuterPosition') - ...
     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
end

%%
set(findall(gcf,'-property','FontSize'),'FontSize',15)
print('-djpeg', '-r0', './Plots/tmp/Prediction.jpg')
print('-dpdf', '-r0', './Plots/tmp/Prediction.pdf')
savefig ./Plots/tmp/Prediction.fig






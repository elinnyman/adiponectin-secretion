
%% Plot prediction
% Plots the prediction, with the values from QuanitifyPrediction"
close all


[~, ind1]=min(values);
[~, ind2]=max(values);
ind=reshape([ind2(1:2:5); ind1(2:2:6)],1,length(ind1));

minVal=min(values(:,2:2:6)); % 'fsk+'; 'fsk-' ; 'noCa+'; 'noCa-';'Iono10+'; 'Iono1.5-'
maxVal=max(values(:,1:2:5));

figure;
set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);

subplot(2,2,3)
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
set(gca, 'XTick', 1:4, 'XTickLabel', {'Control', 'FSK/IBMX', 'FSK/IBMX+BAPTA','FSK/IBMX+Iono'},'XTickLabelRotation',45);

subplot(2,2,4)
errorbar(1:4, [1 3.2 3.21 5.9], [0 0.6 0.5 0.9],'ko','markerface','auto','linewidth', 2)
hold on
title('Experimental Data')
ylabel('Released Adiponectin\newline(fold change)')
set(gca, 'XTick', 1:4, 'XTickLabel', {'Control', 'FSK/IBMX', 'FSK/IBMX+BAPTA','FSK/IBMX+Iono'},'XTickLabelRotation',45);
axis([0.5 4.5 0 10])


param=vp.Parameters(ind,:);

basalValues=exp(param(:,end-2:end));
parameters=basalValues.*change(2:end,:);
colors=reshape([copper(3), copper(3)]',3,6)';
adipoInd=ismember(IQMstates(model),'Adiponectin');

simValues=[];
for k=1:size(param,1)
    
    [simulatedExperiments, allSimData]=SimulateExperiments(param(k,:), 0:1:30*60, design([1; k+1;end]), [basalValues(k,:); parameters(k,:)]);   
    
    simValues=[simValues; allSimData(2).statevalues(2:end,adipoInd)'./allSimData(1).statevalues(2:end,adipoInd)'];
    hold on
end


titles={'FSK/IBMX'; 'BAPTA'; 'Ionomycin'};
for k=1:3
    subplot(2,3,k)
    h=fill([allSimData(1).time(2:end) fliplr(allSimData(1).time(2:end))]/60, [simValues(2*k-1,:) fliplr(simValues(2*k,:))],[255 120 0]./(256+90*(k-1)));
    set(h,'facealpha',1,'EdgeColor','None')
    hold on
    errorbar(30,mean(simValues(2*k-1:2*k,end)), (simValues(2*k-1,end)-simValues(2*k,end))/2,'k','markerface','auto','linewidth',2);
    title(titles{k})
    ylabel('Released Adiponectin\newline(fold change)')
    xlabel('Time (min)')
    axis([0 32 0 9.5])
    
    box off
        set(gca,'FontName', 'Arial')
   
end
set(findall(gcf,'-property','FontSize'),'FontSize',15)
print('-djpeg', '-r0', './Plots/tmp/Prediction.jpg')
print('-dpdf', '-r0', './Plots/tmp/Prediction.pdf')
savefig ./Plots/tmp/Prediction.fig





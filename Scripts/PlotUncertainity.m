function []=PlotUncertainity(validParameters, n, t)
global time
global dgf
global lb
global ub
global wTotalData
global expData
global highCaExpData
global highCaTime
global expData_cAMP
global time_cAMP

if nargin<2
    n=150;
end
if nargin<3
    t=time(end);
end

% 1. Find pararameter sets. First removes all with cost larger than
% chi2-limit, then samples parameters at a linear distance defined by n
validParameters.Cost=validParameters.OriginalCost+(validParameters.CaCost+validParameters.cAMPCost)*wTotalData;

validParameters=unique(validParameters(validParameters.Cost<chi2inv(0.95,dgf-size(validParameters.Parameters,2)),:),'rows');
invalid=find(sum(validParameters.Parameters<repmat(lb,size(validParameters,1),1),2)>0);
invalid=[invalid; find(sum(validParameters.Parameters>repmat(ub,size(validParameters,1),1),2)>0)];
invalid=unique(invalid);
validParameters(invalid,:)=[];

if height(validParameters)>1 % only runs if there exists some valid parameter sets.
    params=validParameters.Parameters;
    selectedParams=[];
    if height(validParameters)>10*n % only samples if there are more parameter sets than 10*n, otherwise uses all parameter sets.
        
        for i=1:size(params,2) % sorts based on parameter i, selects  n parameters with a linear distance.
            params=sortrows(params,i);
            selectedParams=[selectedParams;  params(round(linspace(1,size(params,1),n)),:)];
        end
        validParameters=sortrows(validParameters,2);
        selectedParams=[selectedParams; validParameters.Parameters([1 end],:)];
    else
        selectedParams=validParameters.Parameters;
    end
    params=unique(selectedParams,'rows');
    
    
    % 2. Simulate all selected parameter sets.
    figure(1); % A plot containing all individual simulations
    set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);
    times=0:1:t;
    simulations=table(zeros(size(params,1),length(times)),zeros(size(params,1),length(times)),zeros(size(params,1),length(times)),zeros(size(params,1),length(times)),zeros(size(params,1),length(times)),zeros(size(params,1),length(times)),zeros(size(params,1),length(times)),zeros(size(params,1),length(times))...
        ,'VariableNames',{'Ca_ATP','Ca_noATP','noCa_ATP','noCa_noATP','HiCa_ATP','Ca_ATP_nocAMP','noCa_ATP_nocAMP','noCa_noATP_nocAMP'}); % Define a table to store all simulations in.
    for i = 1:size(params,1)
        param=params(i,:);
        [simulatedExperiments]=SimulateExperiments(param, 0:1:t);
        simulations.Ca_ATP(i,:)=simulatedExperiments.measures(1,:);
        simulations.Ca_noATP(i,:)=simulatedExperiments.measures(2,:);
        simulations.noCa_ATP(i,:)=simulatedExperiments.measures(3,:);
        simulations.noCa_noATP(i,:)= simulatedExperiments.measures(4,:);
        simulations.HiCa_ATP(i,:)= simulatedExperiments.measures(5,:);
        simulations.Ca_ATP_nocAMP(i,:)= simulatedExperiments.measures(6,:);
        simulations.noCa_ATP_nocAMP(i,:)= simulatedExperiments.measures(7,:);
        simulations.noCa_noATP_nocAMP(i,:)= simulatedExperiments.measures(8,:);
        for j=1:6
            subplot(3,2,j)
            hold on
            plot(simulatedExperiments.measures(end,:)./60,simulatedExperiments.measures(j,:)) % Plots all simulations.
        end
    end
    for j=1:6
        subplot(3,2,j)
        hold on
        if j<5
            errorbar(time./60,expData.meanValues(j,:),expData.SEMValues(j,:),'ko','MarkerFaceColor','k')
            axis([0 12.5 -1 24])
        elseif j==5
            errorbar(highCaTime./60,highCaExpData.meanValues,highCaExpData.SEMValues,'ko','MarkerFaceColor','k')
            axis([0 8.5 -6 24])
        elseif j==6
            errorbar(time_cAMP./60,expData_cAMP.meanValues(1,:),expData_cAMP.SEMValues(1,:),'ko','MarkerFaceColor','k')
            axis([0 8.5 -2.5 2.5])
        end
        title(simulatedExperiments.design{j})
        
    end
    
    
    % 3. Plot simulations. This plot is the same as the other figure, but
    % plotted as a filled area instead of individual simulations.
    col=[255 120 0]./256; % orange
    % col=[0 0 256]./256; % blue
    alpha=0.25;
    figure(2);
    hold on
    box off
    
    x=times./60;                  %initialize x array
    subplot(3,2,1)
    y1=max(simulations.Ca_ATP);                      %create first curve
    y2=min(simulations.Ca_ATP);                  %create second curve
    X=[x,fliplr(x)];                %create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %create y values for out and then back
    h1=fill(X,Y,col);                  %plot filled area
    set(h1, 'EdgeColor','None')
    set(gca,'FontName', 'Arial')
    
    subplot(3,2,2)
    y1=max(simulations.Ca_noATP);                      %create first curve
    y2=min(simulations.Ca_noATP);                  %create second curve
    X=[x,fliplr(x)];                %create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %create y values for out and then back
    h2=fill(X,Y,col);                  %plot filled area
    set(h2, 'EdgeColor','None')
    set(gca,'FontName', 'Arial')
    
    subplot(3,2,3)
    y1=max(simulations.noCa_ATP);                      %create first curve
    y2=min(simulations.noCa_ATP);                  %create second curve
    X=[x,fliplr(x)];                %create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %create y values for out and then back
    h3=fill(X,Y,col);                  %plot filled area
    set(h3, 'EdgeColor','None')
    set(gca,'FontName', 'Arial')
    
    subplot(3,2,4)
    y1=max(simulations.noCa_noATP);                      %create first curve
    y2=min(simulations.noCa_noATP);                  %create second curve
    X=[x,fliplr(x)];                %create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %create y values for out and then back
    h4=fill(X,Y,col);                  %plot filled area
    set(h4, 'EdgeColor','None')
    set(gca,'FontName', 'Arial')
    
    subplot(3,2,5)
    y1=max(simulations.HiCa_ATP);                      %create first curve
    y2=min(simulations.HiCa_ATP);                  %create second curve
    X=[x,fliplr(x)];                %create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %create y values for out and then back
    h1=fill(X,Y,col);                  %plot filled area
    set(h1, 'EdgeColor','None')
    set(gca,'FontName', 'Arial')
    
    subplot(3,2,6)
    y1=max(simulations.Ca_ATP_nocAMP);                      %create first curve
    y2=min(simulations.Ca_ATP_nocAMP);                  %create second curve
    X=[x,fliplr(x)];                %create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %create y values for out and then back
    h2=fill(X,Y,col);                  %plot filled area
    set(h2, 'EdgeColor','None')
    set(gca,'FontName', 'Arial')
    
    
    [simulatedExperiments]=SimulateExperiments(validParameters.Parameters(1,:), 0:1:time(end)); % Selects the parameter set with the best agreement with data. 
    validParameters.Cost(1,:)
    set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.6,1]);
    
    for j=1:6 % Plot the simulation with the best agreement with data. 
        subplot(3,2,j)
        hold on
        plot(simulatedExperiments.measures(end,:)./60,simulatedExperiments.measures(j,:),'linewidth',2,'color','k')
        hold on
        if j<5
            errorbar(time./60,expData.meanValues(j,:),expData.SEMValues(j,:),'ko','MarkerFaceColor','k')
            axis([0 12.5 -1 24])
            set(gca,'ytick',[0 10 20], 'xtick',0:2:12)
        elseif j==5
            errorbar(highCaTime./60,highCaExpData.meanValues,highCaExpData.SEMValues,'ko','MarkerFaceColor','k')
            axis([0 8.5 -15 24])
            set(gca,'ytick',[-10 0 10 20], 'xtick',0:2:12)
        elseif j==6
            errorbar(time_cAMP./60,expData_cAMP.meanValues(1,:),expData_cAMP.SEMValues(1,:),'ko','MarkerFaceColor','k')
            axis([0 8.5 -2.5 2.5])
            set(gca,'ytick',[-2 0 2], 'xtick',0:2:12)
        end
        title(simulatedExperiments.design{j})
        
        xlabel('Time (Min)')
        ylabel('\DeltaC/\Deltat (fF/s)')
        box off
        set(gca,'FontName', 'Arial')
    end
    
    cd Plots\tmp
    figure(1)
    set(findall(gcf,'-property','FontSize'),'FontSize',15)
    savefig uncertainty(all).fig
    print -djpeg -r0 uncertainty(all).jpg
    
    figure(2)
    set(findall(gcf,'-property','FontSize'),'FontSize',15)
    savefig uncertainty.fig
    print -dpdf -r0 uncertainty.pdf
    print -djpeg -r0 uncertainty.jpg
    
    cd ..\..
end
end



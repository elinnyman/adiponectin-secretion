function [  ] = PlotInterval( validParameters, n)
global dgf
global lb
global ub
global wTotalData

%% Calculate the endocytosis contribution to the recorded signal. 
endo=validParameters.ExoEndo(:,1:2:5);
exo=validParameters.ExoEndo(:,2:2:6);
validParameters.ExoEndo=endo./(endo+exo);

if nargin==2
    validParameters.ExoEndo=round(validParameters.ExoEndo.*10^n)./10^n; % Rounds to the nth decimal. 
end

%% Removes invalid parameter sets. 
validParameters.Cost=validParameters.OriginalCost+(validParameters.CaCost+validParameters.cAMPCost)*wTotalData; % Updates the total cost to be the sum of all "subcosts" (original + high Ca + 0 cAMP).
invalid=find(sum(validParameters.Parameters<repmat(lb,size(validParameters,1),1),2)>0);
invalid=[invalid; find(sum(validParameters.Parameters>repmat(ub,size(validParameters,1),1),2)>0)];
invalid=unique(invalid);
validParameters(invalid,:)=[];
validParameters=sortrows(validParameters,'Cost');

%% Generate the profile-likelihood plots, and the interval figure. 
titles={'Ca+ATP','Ca+noATP','10xCa+ATP'};
if ~isempty(validParameters)
    for i=1:size(validParameters.ExoEndo,2)
        figure(i)
        hold on
        set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,1,1]);
        params=validParameters.ExoEndo(:,i).*100;
        [~,ind]=unique(params);
        plot(params(ind),validParameters.Cost(ind),'.','linewidth',3,'color',[255 120 0]./256)
        xlabel('Endocytosis contribution')
        ylabel('-error')
        hold on
        axis([0 100 20 36])
        
        limit=chi2inv(0.95,dgf);
        plot([min(params) max(params)],[limit limit])
        
        limit=chi2inv(0.95,dgf-size(validParameters.Parameters,2));
        plot([0 100],[limit limit],'k','linewidth',2)
        set(gca,'ytick',[20 25 30 35])
        figure(i)
        title(titles{i})
        box off
        set(gca,'FontName', 'Arial')
        set(findall(gcf,'-property','FontSize'),'FontSize',15)
       
        name=['PPL' num2str(i)];
        savefig(['./Plots/tmp/' name '.fig'])
        print('-dpdf', '-r0', ['./Plots/tmp/' name '.pdf'])
        print('-djpeg', '-r0', ['./Plots/tmp/' name '.jpg'])
        
        
        optVal=round(min(validParameters.ExoEndo(validParameters.Cost==min(validParameters.Cost),i))*10000)/100;
        minVal=round(min(validParameters.ExoEndo(validParameters.Cost<limit,i))*10000)/100;
        maxVal=round(max(validParameters.ExoEndo(validParameters.Cost<limit,i))*10000)/100;
        fprintf('%i: %.2f %%  (%.2f - %.2f / %.2f +- %.2f %%)  \n', i,  optVal, minVal, maxVal, (maxVal+minVal)/2, (maxVal+minVal)/2-minVal)
        
        figure(4)
        errorbar(i,mean([minVal maxVal]), (maxVal-minVal)/2,'color',[255 120 0]./256,'linewidth',3)
        hold on;
        
    end
    figure(4)
    set(gca, 'XTick', 1:3, 'XTickLabel', titles,'XTickLabelRotation',45,'YTick',0:20:100);
    axis([0.5 3.5 0 100]);
    ylabel('% endocytosis contribution')
    box off
    set(gca,'FontName', 'Arial')
    set(findall(gcf,'-property','FontSize'),'FontSize',15)
    set (gcf, 'PaperPositionMode','auto','Units', 'normalized', 'outerposition', [0,0,0.35,0.6]);
    print('./Plots/tmp/Interval.pdf','-dpdf')
    
    
    disp(['scale: ' num2str(min(validParameters.scale(validParameters.Cost<limit))) ' - ' num2str(max(validParameters.scale(validParameters.Cost<limit)))]);
else
    disp('no valid params')
end
end


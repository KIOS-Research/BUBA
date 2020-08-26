function [ h ] = PaperCLBoundsAndMCS(d,dispname,stateEstim,sensorNodesCLind,saveResults)

% function [ h ] = PlotCLBoundsAndMCS(CLmin, CLmax, CLmean, cnmin, cnmax)
%Description: Plots first set of bounds (internally, filled) with second 
%             set of bounds (outwards,lines) 

%% Assign variables and correct timeseries
cnmin = stateEstim.NodeQualityLB(:,sensorNodesCLind)';
cnmax = stateEstim.NodeQualityUB(:,sensorNodesCLind)';
CLminShort = stateEstim.NodeQualityLBmcs(:,sensorNodesCLind)';
CLmaxShort = stateEstim.NodeQualityUBmcs(:,sensorNodesCLind)';
CLmean=[];
n=(size(cnmin,2)-1)/(size(CLminShort,2)-1);
CLmin=[kron(CLminShort(:,1:end-1),ones(1,n)) CLminShort(:,end)];
CLmax=[kron(CLmaxShort(:,1:end-1),ones(1,n)) CLmaxShort(:,end)];
    

%% Check inputs
if any([isempty(CLmin) isempty(CLmax) isempty(cnmin) isempty(cnmax)])
    disp(sprintf('\nerror: Inputs are empty.'));
    return
end

%% create x-axis time table
time = linspace(0,stateEstim.Time(end)/3600,size(cnmin,2));
if length(time)==1; Colr ='k*'; else Colr ='k'; end

%% Chlorine bounds CYDMA
fontsize=12;
if contains(dispname,'Hanoi')
    maxplotperfig=3;
    maxcolumns = 3;
    yaxislimits = [0.22, 0.4];%[0.33, 0.4];
    figuresize = [0 0 0.95 0.5];
    legendPos='SouthWest';
    xlabelminInd=0;
    if contains(dispname,'wall')
        pdfName='results/Hanoi_wall';
        figName='results/Hanoi_wall.fig';
    else
        pdfName='results/Hanoi';
        figName='results/Hanoi.fig';
    end
elseif contains(dispname,'CY_DMA')
    maxplotperfig=9;
    maxcolumns = 3;
    yaxislimits = [0, 0.4];%[0.33,0.4];
    figuresize = [0 0 0.95 1];
    legendPos='SouthEast';
    xlabelminInd=7;
    if contains(dispname,'wall')
        pdfName='results/CYDMA_wall';
        figName='results/CYDMA_wall.fig';
    else
        pdfName='results/CYDMA';
        figName='results/CYDMA.fig';
    end
end

nodenum=length(sensorNodesCLind);
plotperfig = min(maxplotperfig,nodenum);
columns = min(maxcolumns,nodenum);
lines = ceil(plotperfig/columns);
for j=1:ceil(nodenum/plotperfig)
    f = figure('units','normalized','outerposition',figuresize);
    for i=1:plotperfig
        node=i+(j-1)*plotperfig;
        if node>nodenum; break; end
        subplot(lines,columns,i);
        h(1)=plot(time,cnmin(node,:),Colr,'linewidth',1.5);
        set(gca,'fontsize',fontsize);
        hold all
        plot(time,cnmax(node,:),Colr,'linewidth',1.5)
        if ~isempty(CLmean)
            plot(time,CLmean(node,:))
        end
        h(2)=jbfill(time,CLmin(node,:),CLmax(node,:),'r','r');
        nodeID=d.getNodeNameID(sensorNodesCLind(node));
        title(['Node ',nodeID{1}])
        if i>=xlabelminInd
        xlabel('Time (hours)')
        end
        ylabel('CL (mg/L)')
        axis tight
        if i==1
            legend(h([1 2]), 'BUBA bounds','MCS bounds','Location',legendPos)
            
        end
        grid on
        if ~isempty(yaxislimits)
            set(gca,'ylim',yaxislimits)
        end
    end
    tightfig;
end
if saveResults
print(pdfName,'-dpdf')
savefig(f,figName)
end

%% Threshold violations
% nodenum=size(cnmin,1);
% plotperfig=12;
% for j=1:ceil(nodenum/plotperfig)
% figure('units','normalized','outerposition',[0 0 0.95 1])
% % figure('units','normalized','outerposition',[0 0 0.6 0.6])
% for i=1:plotperfig
%     node=i+(j-1)*plotperfig;
%     if node>nodenum; break; end
%     subplot(3,4,i)
%     minViol = CLmin(node,:)-cnmin(node,:);
%     minViol(minViol>0)=0;
%     maxViol = cnmax(node,:)-CLmax(node,:);
%     maxViol(maxViol>0)=0;
%     viol = abs(minViol + maxViol);
%     plot(time,viol,Colr,'linewidth',1.5)
%     set(gca,'fontsize',fontsize)
%     title(['Violations at node ',d.getNodeNameID(sensorNodesCLind(node))])%num2str(node)])
%     xlabel('Time (quality steps)')
%     ylabel('Violations (mg/L)')
%     axis tight
%     grid on
% end
% % keyboard
% end


end
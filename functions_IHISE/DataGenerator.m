function [ nodeTimeSeries, Qepa, Hepa, LinkStatus ] = DataGenerator( d )
%DEMANDGENERATOR 

%% Epanet Simulation

%calculate
EPAall=d.getComputedTimeSeries;
Hepa=EPAall.Head;
Qepa=EPAall.Flow;
LinkStatus=EPAall.Status;
nodeTimeSeries=EPAall.Demand;

%process
resIndex=d.NodeReservoirIndex;
tankIndex=d.NodeTankIndex;
nodeTimeSeries(:,[resIndex, tankIndex])=Hepa(:,[resIndex, tankIndex]);
nodeTimeSeries=nodeTimeSeries';
Hepa=Hepa';
Qepa=Qepa';
LinkStatus=LinkStatus';

% Ktemp = (EPAall.HeadLoss./(abs(EPAall.Flow).^2));
% Kepa = Ktemp(5,:);
% Kepa = Kepa';
end


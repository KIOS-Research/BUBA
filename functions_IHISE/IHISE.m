function [Qlower,Qupper,hlower,hupper, info] = IHISE(inpname,Kunc,qunc,Qmeas,Hmeas)
%% IHISE

%% Load Network
d=epanet(inpname);

%% Insert uncertainty
% Kunc=0.02; %pipe parameter uncertainty
% qunc=0.02; %demand uncertainty

%% Define simulation times
simTime=d.getTimeSimulationDuration;
hydStep=d.getTimeHydraulicStep;
if simTime==0; simTime=hydStep;end %In case simTime is zero, do one simulation
times=struct('simTime',simTime,'hydStep',hydStep);

%% Generate node time series and simulate with epanet
% nodeTimeseries (n_n x simSteps) contains all the nodes sorted by index
% demand nodes have the node demand for every simulation step
% reservoir nodes have the total head for every simulation step
% tank nodes have the  the total head (elevation + tanklevel) for every simulation step
[ nodeTimeSeries, Qepa, Hepa, LinkStatus ] = DataGenerator( d );

%check for Epanet negative pressures
elev=repmat(double(d.getNodeElevations'),[1,size(Hepa,2)]);
if any(any((Hepa-elev)<0))
    error('Epanet returned negative pressures')
end

%% Run IHISE
info = struct('iterPerTimeStep',[],'execTimePerTimeStep',[],'execTimeTotal',[]);
startIHISE = tic;
[ Qlower, Qupper, hlower, hupper, info] = IHISE_TimeSteps(d, nodeTimeSeries, times, Kunc, qunc, Qepa, Hepa, LinkStatus, Qmeas,  Hmeas, info );
timeIHISE = toc(startIHISE);
info.execTimeTotal = timeIHISE;

%%
d.unload; %unloads libraries and deletes temp files

end

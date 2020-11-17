try 
d.unload
catch ERR
end 
fclose all;clear class;clear;clc;close all;
addpath(genpath(pwd));
% load_paths()

%% Simulation settings:
ihiseOverwrite = 0; % requires Gurobi and to uncomment load_paths()
mcsOverwrite = 0;
bubaOverwrite = 1;
allOverwrite = 1;
mcsFlowBounds = 1;
allSensorNodes = 0; 
mcs_iter = 30000; 
displayResults = 1;

%% Choose Network:
[inpname,dispname] = enterNetwork([]);

%% Network uncertainty
qunc=0.1; %demand uncertainty
Kunc=0.1; %pipe parameter uncertainty
Kdunc=0.1; %decay rate uncertainty

%% Define chlorine sensor locations:
d=epanet(inpname);
if ~allSensorNodes
    if contains(dispname,'Hanoi')
        sensorNodesCLid = {'5','10','24'}; %quality sensor ids
    elseif contains(dispname,'CY_DMA')
        sensorNodesCLid = {'3','17','28','41','44','55','65','82','90'}; %quality sensor ids
    else
        error('Unknown sensor locations in given network.')
    end
else
    sensorNodesCLid = d.getNodeNameID;
end

%% Simulation times:
th = double(d.getTimeHydraulicStep); 
tq = double(d.getTimeQualityStep);
tstep = tq;

%% State estimation:
stateEstim = d.getComputedTimeSeries;
d.unload

%%% get quality state data in quality time step:
if (th~=tq)
    d=epanet(inpname);
    d.setTimeReportingStep(tq);
    allParameters=d.getComputedTimeSeries;
    stateEstim.NodeQualityTQ = allParameters.NodeQuality;
    clearvars allParameters
    d.unload
end

%% IHISE simulation:
%%% check if file exists:
filenameIHISE = [pwd,'\simulations\IHISE_',dispname,'.mat'];
if isfile(filenameIHISE) && ~ihiseOverwrite
    % File exists, load:
    load(filenameIHISE)
    disp('Saved IHISE simulation loaded!')
elseif ihiseOverwrite
    % File does not exist, calculate:
    [Qlower, Qupper, hlower, hupper, IHISEinfo] = IHISE(inpname,Kunc,qunc,[],[]);%Generates and saves bounds
    save(filenameIHISE,'Qlower', 'Qupper', 'hlower', 'hupper', 'IHISEinfo', 'Kunc', 'qunc')
else
    % requested to load but file does not exist:
    error('IHISE file does not exist!')
end

%%% Include hydraulic state bounds in state estimation:
stateEstim.FlowUB = Qupper';
stateEstim.FlowLB = Qlower';
stateEstim.HeadUB = hupper';
stateEstim.HeadLB = hlower';
clearvars Qupper Qlower hupper hlower

%% Calculate or load MCS bounds
%%% check if file exists:
filenameMCS = [pwd,'\simulations\MCS_',dispname,'.mat'];
if isfile(filenameMCS) && ~mcsOverwrite
    % File exists, load:
    load(filenameMCS)
    disp('Saved MCS loaded!')
elseif mcsOverwrite
    % File does not exist, calculate:
    [Qmin, Qmax, Qmean, hmin, hmax, hmean, CLmin, CLmax, CLmean, MCSinfo] = MCSCL2(mcs_iter,tstep,inpname,Kunc,qunc,Kdunc);
    save(filenameMCS,'Qmin', 'Qmax', 'Qmean', 'hmin', 'hmax', 'hmean', 'CLmin', 'CLmax', 'CLmean', 'MCSinfo')
else
    % requested to load but file does not exist:
    error('MCS file does not exist!')
end

%%% Include MCS state bounds in state estimation:
stateEstim.FlowUBmcs = Qmax';
stateEstim.FlowLBmcs = Qmin';
stateEstim.HeadUBmcs = hmax';
stateEstim.HeadLBmcs = hmin';
stateEstim.NodeQualityUBmcs = CLmax';
stateEstim.NodeQualityLBmcs = CLmin';
clearvars Qmax Qmin hmax hmin CLmin CLmax Qmean hmean CLmean

%% BUBA simulation:
%%% check if file exists:
filenameBUBA = [pwd,'\simulations\BUBA_',dispname,'.mat'];

if isfile(filenameBUBA) && ~bubaOverwrite
    % File exists, load:
    load(filenameBUBA)
    disp('Saved BUBA simulation loaded!')
elseif bubaOverwrite
    % Calculate BUBA:
    [cnmin, cnmax, sensorNodesCLind] = BUBA(inpname,sensorNodesCLid,stateEstim,Kdunc,mcsFlowBounds); % non-varying decay rate
    save(filenameBUBA,'cnmin', 'cnmax','sensorNodesCLind', 'Kunc', 'qunc', 'Kdunc','mcsFlowBounds')
else
    % requested to load but file does not exist:
    error('BUBA file does not exist!')
end

%%% Include quality state bounds in state estimation:
stateEstim.NodeQualityUB = cnmax';
stateEstim.NodeQualityLB = cnmin';
clearvars cnmin cnmax

%% Save all and unload
if allOverwrite
    filename = [pwd,'\simulations\ALL_',dispname,'.mat'];
    d=epanet(inpname);
    save(filename)
    d.unload
end

%% Display results
if displayResults
    run DisplayResults.m
end

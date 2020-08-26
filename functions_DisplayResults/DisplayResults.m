try 
d.unload
catch ERR
end 
fclose all;clear class;clear;clc;close all;
addpath(genpath(pwd));
% load_paths()

%% Choose Network
% inpname = enterNetwork([]);
dataname = enterData([]);
load([pwd,dataname])
d=epanet(inpname);

%% Paper resutls
%%% Define parameters:
allSensorNodes=0;
saveResults=1;
%%% Define sensor nodes:
if ~allSensorNodes
    if contains(dispname,'Hanoi')
        sensorNodesCLid = {'5','10','24'};
        paperSensorNodesCLind = d.getNodeIndex(sensorNodesCLid); %quality sensor ids
    elseif contains(dispname,'CY_DMA')
        sensorNodesCLid = {'3','17','28','41','44','55','65','82','90'}; %quality sensor ids
        paperSensorNodesCLind = d.getNodeIndex(sensorNodesCLid);
    else
        error('undefined network sensors')
    end
else
    paperSensorNodesCLind = sensorNodesCLind;
end
%% Plot:
PaperCLBoundsAndMCS(d,dispname,stateEstim,paperSensorNodesCLind,saveResults)

%%
d.unload
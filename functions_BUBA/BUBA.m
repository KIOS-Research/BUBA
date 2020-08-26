function [cnmin, cnmax, sensorNodesCL] = BUBA(inpname,sensorNodesCLid,stateEstim,Kdunc,mcsFlowBounds)

%% Network quality setting check:
NetworkQualitySettingsCheck(inpname);

%%
bita=Kdunc; 
Q = stateEstim.Flow;
if mcsFlowBounds
    Ql = stateEstim.FlowLBmcs;
    Qu = stateEstim.FlowUBmcs;
else
    Ql = stateEstim.FlowLB;
    Qu = stateEstim.FlowUB;
end
H = stateEstim.Head;
CL = stateEstim.NodeQualityTQ;

d=epanet(inpname);
sim_time=double(d.getTimeSimulationDuration);
sensorNodesCL = double(d.getNodeIndex(sensorNodesCLid)); %sensor nodes

%% Pipe Area Ar, Pipe length L
diam=d.getLinkDiameter/1000; %in meters
radius = diam/2;
Ar = pi*(diam./2).^2; %pipe area= pi*r^2 in (meters^2)
L = d.getLinkLength; %pipe lenght in meters

%% Tank and reservoir parameters
TankID=d.getNodeTankIndex;
ReservoirID=d.getNodeReservoirIndex;
TankV=(H(:,TankID)- ones(size(H(:,1)))*d.NodeElevations(TankID)).*(pi*ones(size(H(:,1)))*(d.NodeTankDiameter(TankID)/2).^2); %

%% Create incidence matrix (topology graph)
node_node_link_index=[d.NodesConnectingLinksIndex d.LinkIndex'];
node_node_link_name=[d.NodesConnectingLinksID d.LinkNameID'];
Ain=zeros(d.NodeCount,d.LinkCount); %size nodes x pipes
for i=1:length(d.NodesConnectingLinksIndex)
    Ain(node_node_link_index(i,1),node_node_link_index(i,3))=1;
    Ain(node_node_link_index(i,2),node_node_link_index(i,3))=-1;
end

%% Kb, Kw, K, Ke, Ku, Kl - Reaction rates
%%%% Kb - Bulk reaction coefficient in hours
Kb = d.getLinkBulkReactionCoeff./24; % Bulk reaction coefficient (1/hour) 
d.QualityReactionCoeffBulkUnits; % Bulk reaction coefficient units

%%%% Kw - Wall reaction coefficient in hours
Kw = d.getLinkWallReactionCoeff./24; %Wall reaction coefficient (m/hour)
d.QualityReactionCoeffWallUnits; % Wall reaction coefficient units

%%%% Ktank - Tank reaction coefficient in hours
KTank=d.getNodeTankBulkReactionCoeff;
KTank = KTank./24;
KlTank=KTank - bita*KTank;
KuTank=KTank + bita*KTank;

%%%% Kf - Calculation of mass tranfer coefficient
Kfw = (4./diam).*Kw; % page 44 epanet manual

%%% K - Total chlorine decay rate:
K = (Kb + Kfw);

%%%% Kl, Ku - Bounds on total chlorine decay rate 
Kl = K + bita.*K; %negative plus negative = most negative
Ku = K - bita.*K; %negative minus negative = less negative

%% Set simulation steps and Convert simulation time steps from seconds to hours 
th = double(d.getTimeHydraulicStep)/3600; %Hydraulic time step in hours
tq = double(d.getTimeQualityStep) /3600; %Quality Step in hours
start_step=1;
end_step=(sim_time/3600/tq)+1; %Calculate last time step

%% Replicate missing hydraulic steps
% if th~=tq
if size(CL,1)~=size(Ql,1)
n=th/tq;
TankV=[repelem(TankV(1:end-1,:)',1,n) TankV(end,:)']';
n=(end_step-1)/(size(Ql,1)-1);
Ql=[repelem(Ql(1:end-1,:)',1,n) Ql(end,:)']';
Qu=[repelem(Qu(1:end-1,:)',1,n) Qu(end,:)']';
end

%% Initialize calculated bounds and set as global variable
%initialize all node quality:
global c_calc
c_calc=NaN(d.NodeCount,2,end_step); 

%set initial quality of nodes:
c_calc(:,:,1)=[CL(1,:)' CL(1,:)'];
c_calc(:,:,2)=[CL(2,:)' CL(2,:)'];

%find input nodes:
inputN=[];
type=d.getNodeSourceType;
for i=1:length(type)
    if strcmp(type(i),'SETPOINT')
        inputN = [inputN i];
    end
end

%set the quality of input nodes:
for i = inputN
c_calc(i,1,1:end_step)=CL(1:end_step,i); %input node minimum
c_calc(i,2,1:end_step)=CL(1:end_step,i); %input node maximum
end

%% Constant Parameter and Calculated hydraulic states struct:
times = struct('QualityStep',tq, 'HydraulicStep',th, 'SimulationTime',sim_time);
links = struct('FlowUpper',Qu, 'FlowLower',Ql, 'Length',L, 'Area',Ar, 'DecayRateUpper',Ku, 'DecayRateLower',Kl);
tanks = struct('ID',TankID, 'Volume',TankV, 'DecayRateUpper',KuTank, 'DecayRateLower',KlTank);
nodes = struct('TankInfo',tanks, 'ReservoirID',ReservoirID, 'ActualClConcen',CL);
sim_param = struct('time',times, 'links',links, 'tanks',tanks, 'nodes',nodes, 'IncidenceMat',Ain);

%% Find output node min-max concentration
i=1;
starttime=tic;
Max_Loop_Time=0;
cnmin=zeros(1,end_step);
cnmax=zeros(1,end_step);
nn=d.getNodeCount;
sensor_num_all = length(sensorNodesCL);
sensor_num = 1;
for Na=sensorNodesCL
for kc=start_step:end_step
    loopstarttime=tic;
    [cnmin(kc), cnmax(kc)] = Algorithm_5(Na,kc,sim_param);
    clc
    %%%%%%%%%%%%%%%%%% Print times  %%%%%%%%%%%%%%
    Elapsed_Time=toc(starttime);
    loop_Elapsed_Time(i)=toc(loopstarttime);
    if i>10
    Max_Loop_Time=max(Max_Loop_Time,loop_Elapsed_Time(i));
    end
    fprintf('The program elapsed time is:  %i:%.2i  minutes:seconds \n\n', (floor(Elapsed_Time/60)),floor(mod(Elapsed_Time,60)));
    fprintf('The current computation step time is:  %.2f  seconds\n\n', loop_Elapsed_Time(i));
    Average_Loop_Time=mean(loop_Elapsed_Time);
    fprintf('The average computation step time is:  %.2f  seconds\n\n', Average_Loop_Time);
    fprintf('The maximum computation step time is:  %.2f  seconds\n\n', Max_Loop_Time);
    fprintf('Sensor number: --> %i <-- out of  %i  sensors\n\n', sensor_num,sensor_num_all);
    fprintf('Time step is : --> %i <-- out of  %i  steps\n\n', kc,end_step);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
i=i+1;
end
sensor_num = sensor_num+1;
end
fprintf('\nBUBA Simulation Finished\n\n');

%% reshape results:
cnmax=reshape(c_calc(:,2,:),size(c_calc,1),size(c_calc,3));
cnmin=reshape(c_calc(:,1,:),size(c_calc,1),size(c_calc,3));

%%
d.unload
end


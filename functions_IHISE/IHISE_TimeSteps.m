function [ Qlower, Qupper, hlower, hupper, info] = IHISE_TimeSteps(d, nodeTimeSeries, times, Kunc, qunc, Qepa, Hepa, LinkStatus, Qmeas, Hmeas, info )

%% Load neccesary classes
intv=intvclass2();

%% Set times
simTime=times.simTime; %simulation time in seconds
% simSteps=floor(simTime/times.hydStep);
simSteps=size(nodeTimeSeries,2);
restankIndex=[d.getNodeTankIndex d.getNodeReservoirIndex];
restankIndex(find(restankIndex==0))=[];

%% Formulate demand matrix and reservoir and tank level matrix for IHISE
%nodeTimeSeries: (n x simsteps), lines in order with node index

%% Pipe resistance coefficients (K, n) and uncertainty:
[ K, n ] = Pipe_Coefficients( d );
%check pipe resistance coefficients
Ktest=K;
Ktest(d.LinkPumpIndex)=[];
if any(Ktest<=0)
    error('Invalid pipe resistance coefficient')
end
clear Ktest
%no pumps
%Random bounds within Kunc%
Rand=2*rand(length(K),1);
Kl=K-Rand*Kunc.*K;
Ku=K+(2-Rand)*Kunc.*K;
Kbnd=[Kl Ku];

%Equal bounds within Kunc%
Kl=K-Kunc*K;
Ku=K+Kunc*K;
Kbnd=[Kl Ku];

%% Initialize bounds
Qlower=zeros(d.LinkCount,simSteps);
Qupper=zeros(d.LinkCount,simSteps);
hlower=zeros(d.NodeCount,simSteps);
hlower(restankIndex,:)=[];
hupper=zeros(d.NodeCount,simSteps);
hupper(restankIndex,:)=[];

%% Pump coefficients
%Calculate pump curve function coefficients
[ Pcoef ] = Pump_Coefficients( d );

%% IHISE Extended Time Simulation 
starttime=tic;
for simStep=1:simSteps

    loopstarttime=tic;
    
    %% External demands and uncertainty:
    qext=nodeTimeSeries(:,simStep);
    qext(restankIndex)=[];
    
%     %Random bounds within qunc%
%     Rand=2*rand(length(qext),1);
%     ql=qext-Rand*qunc.*abs(qext);
%     qu=qext+(2-Rand)*qunc.*abs(qext);
%     qbnd=[ql, qu];
    
    %Enallax bounds panw katw
    %To see if bounds tighten due to making most solutions infeasible
%     Enalbin1=ones(size(qext));
%     Enalbin1(find(mod(find(Enalbin1==1),2)==0))=0;
%     ql=qext-Enalbin1*2*qunc.*qext;
%     Enalbin2=zeros(size(qext));
%     Enalbin2(find(mod(find(Enalbin2==0),2)==0))=1;
%     qu=qext+Enalbin2*2*qunc.*qext;
%     qbnd=[ql, qu];
    
    %Equal bounds within qunc
    ql=qext-qunc*qext;
    qu=qext+qunc*qext;
    qbnd=[ql, qu];
    
    %% Reservoir and tank levels:
    %restankLevels is a vector size n_n of zeros with values only at
    %reservoir and tank index
    restankLevels = nodeTimeSeries(:,simStep);
    temp=ones(size(restankLevels));
    temp(restankIndex)=0;
    restankLevels(temp==1)=0;
    linkStatusStep = LinkStatus(:,simStep);
    
    %% Conventional State Estimation
%     [Q, h] = Snapshot_Hyd_SE(d, restankLevels, linkStatusStep, K, n, Pcoef, qext);  
%     h=double(h);
%     % check if snapshot solution agrees with Epanet solution
%     epsilon=2;
%     if any(abs(Qepa(:,simStep)-Q)>epsilon) || any(abs(Hepa(1:length(h),simStep)-h)>epsilon)
%         [Qepa(:,simStep) Q];
%         [Hepa(1:length(h),simStep) h];
%         maxFlowDeviation = max(abs(Qepa(:,simStep)-Q))
%         maxHeadDeviation = max(abs(Hepa(1:length(h),simStep)-h))
%         warning('Newton method does not agree with EPANET.')
%         pause(2)
%         Q=Qepa(:,1);h=Hepa(1:size(h,1),1);%replace with EPANET estimates
%     end
        
    Q=Qepa(:,1);h=Hepa(1:size(ql,1),1);%replace with EPANET estimates
    
    %% Find initial bounds
    [ Qbnd, hbnd ] = Initial_Bounds( d, intv, n, restankLevels, linkStatusStep, Kbnd, qbnd, Pcoef, Q, h );
    
    %%Extra measurements:
    QSenUnc = 0.02;
    for i=Qmeas
        Qbnd(i,1)=Q(i) - QSenUnc*abs(Q(i));
        Qbnd(i,2)=Q(i) + QSenUnc*abs(Q(i));
    end
    HSenUnc = 0.002;
    elev = double(d.getNodeElevations); %remove elevation offset from sensor uncertaity
    elev = elev(1:length(h))';
    for i=Hmeas
        hbnd(i,1)=h(i) - HSenUnc*(h(i)-elev(i));
        hbnd(i,2)=h(i) + HSenUnc*(h(i)-elev(i));
    end
    
    % check if the feasible solution is within the initial bounds
    if (any(Qbnd(:,1)>Q) || any(Qbnd(:,2)<Q) || any(hbnd(:,1)>h) || any(hbnd(:,2)<h))
        [Qbnd(:,1) Q Qbnd(:,2)];
        [hbnd(:,1) h hbnd(:,2)];
        warning('Known feasible solution is not within the initial bounds.')
        pause(2)
        % fix wrong bounds using estimates
        Qbnd((find(Qbnd(:,1)>Q)),1) = Q(find(Qbnd(:,1)>Q))-100*abs(Q(find(Qbnd(:,1)>Q)));
        Qbnd((find(Qbnd(:,2)<Q)),2) = Q(find(Qbnd(:,2)<Q))+100*abs(Q(find(Qbnd(:,2)<Q)));
        hbnd((find(hbnd(:,1)>h)),1) = h(find(hbnd(:,1)>h))-10*abs(h(find(hbnd(:,1)>h)));
        hbnd((find(hbnd(:,2)<h)),2) = h(find(hbnd(:,2)<h))+10*abs(h(find(hbnd(:,2)<h)));        
    end
    
    %% Iterations to find bounds on Q and h
    startIHISE=tic;
    
    [ Q, h, info ] = IHISE_Iterations( d, n, restankLevels, linkStatusStep, Kbnd, qbnd, Qbnd, hbnd, Pcoef, info);
    
    execTime = toc(startIHISE);
    info.execTimePerTimeStep = [info.execTimePerTimeStep execTime];
    %%%%%%%%%%%%%%%%%% Print times  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc
    Elapsed_Time=toc(starttime);
    loop_Elapsed_Time(simStep)=toc(loopstarttime);
    fprintf('*** Iterative Hydraulic Interval State Estimation (IHISE) *** \n\n');
    fprintf('The program elapsed time is:  %i:%.2i  minutes:seconds \n\n', (floor(Elapsed_Time/60)),floor(mod(Elapsed_Time,60)));
    fprintf('The current  loop   time is:  %.2f  seconds\n\n', loop_Elapsed_Time(simStep));
    Average_Loop_Time=mean(loop_Elapsed_Time);
    fprintf('The average  loop   time is:  %.2f  seconds\n\n', Average_Loop_Time);
    fprintf('Time step is: --> %i <-- of   %i   steps\n\n', simStep,simSteps);
    est_rem_time=Average_Loop_Time*(simSteps-simStep);
    est_rem_time_min=floor(est_rem_time/60);
    est_rem_time_sec=floor(mod(est_rem_time,60));
    fprintf('The estimated remaining time is:  %i:%.2i  minutes:seconds \n\n', est_rem_time_min,est_rem_time_sec);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Save bounds timeseries
    Qlower(:,simStep)=Q(:,1);
    Qupper(:,simStep)=Q(:,2);
    hlower(:,simStep)=h(:,1);
    hupper(:,simStep)=h(:,2);

end


end
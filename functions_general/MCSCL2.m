function [Qmin, Qmax, Qmean, hmin, hmax, hmean, CLmin, CLmax, CLmean, info] = MCSCL2(max_iter,tstep,inpname,Kunc,qunc,Kdunc)
%%Monte-Carlo Simulations for hydraulics and quality
    
%% Fresh start:
try 
d.unload
catch ERR
end 

%% Load Network
d=epanet(inpname);
strstart = max(strfind(inpname,'\'))+1;
dispname=inpname(strstart:end-4);

%% Set times
d.setTimeReportingStep(tstep);

%% Demand Bounds
temp = d.NodeBaseDemands;
qext=temp{1};
dem_length=length(qext);
ql=qext-qunc*qext;
qu=qext+qunc*qext;

%% Rougness coefficient bounds
L=d.getLinkLength';
L_length = length(L);
Ll=L-Kunc*L;
Lu=L+Kunc*L;

%% Bulk reaction coefficients bounds
KTank=d.getNodeTankBulkReactionCoeff;
KlTank=KTank - Kdunc*KTank;
KuTank=KTank + Kdunc*KTank;
K = d.getLinkBulkReactionCoeff;
K_length = length(K);
Kl = K + Kdunc.*K; %negative plus negative = most negative
Ku = K - Kdunc.*K; %negative minus negative = less negative

%% Wall reaction coefficients bounds
Kw = d.getLinkWallReactionCoeff;
Kwl = Kw + Kdunc.*Kw; %negative plus negative = most negative
Kwu = Kw - Kdunc.*Kw; %negative minus negative = less negative

%% Inititalize:
Sim = d.getComputedTimeSeries_ENepanet;
Qmin = Sim.Flow';
Qmax = Sim.Flow';
Qmean = Sim.Flow';
hmin = Sim.Head';
hmax = Sim.Head';
hmean = Sim.Head';
CLmin = Sim.NodeQuality';
CLmax = Sim.NodeQuality';
CLmean = Sim.NodeQuality';

%% Monte-Carlo Solve Hydraulics with EPANET
startTime = tic;
iter = 1;
while iter<max_iter  
    iter = iter+1;
    disp(['Monte-Carlo simulation : ',num2str(iter)])
    
    %Vary demands
    qext=ql+rand(1,dem_length).*(qu-ql);
    d.setNodeBaseDemands({qext});
    
    %Vary model uncertainty (pipe lengths)
    if Kunc>0.0001
    L=Ll+rand(L_length,1).*(Lu-Ll);
    d.setLinkLength(L)
    end
    
    %Vary model uncertainty (bulk reaction)
    K = Kl+rand(1,K_length).*(Ku-Kl);
    d.setLinkBulkReactionCoeff(K);
    
    %Vary model uncertainty (wall reaction)
    Kw = Kwl+rand(1,K_length).*(Kwu-Kwl);
    d.setLinkWallReactionCoeff(Kw);
    
    %%Simulate
    Sim = d.getComputedTimeSeries_ENepanet;
    Q = Sim.Flow';
    h = Sim.Head';
    CL = Sim.NodeQuality';
    
    %%Save min, mean and max state:
    Qmin=min(Qmin,Q);
    Qmax=max(Qmax,Q);
    Qmean = ((iter-1)*Qmean + Q)/iter;
    hmin=min(hmin,h);
    hmax=max(hmax,h);
    hmean = ((iter-1)*hmean + h)/iter;
    CLmin=min(CLmin,CL);
    CLmax=max(CLmax,CL);
    CLmean = ((iter-1)*CLmean + CL)/iter; 
    
end
if exist('startTime','var')
    timeMCS = toc(startTime);
else
    timeMCS = NaN;   
end
info = struct('MCStime',timeMCS,'MCSiter',max_iter);
          
%% Unload
d.unload
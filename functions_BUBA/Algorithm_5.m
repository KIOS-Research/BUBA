function [ cnmin, cnmax ] = Algorithm_5(Na, kc, sim_param)
%ALGORITHM_5 
% Finds chlorine concentration bounds of a general node in a water network

% Finds the min-max concentration cnmin-cnmax of node Na at time instant kc.

%% Extract function parameters:
global c_calc

%% if values have already been calculated:
if (~isnan(c_calc(Na,1,kc)) && ~isnan(c_calc(Na,2,kc)))
    cnmin = c_calc(Na,1,kc);
    cnmax = c_calc(Na,2,kc);
    return
end

%% if time step exceeds initial time step:
if kc<=1
    cnmin=c_calc(Na,1,1);
    cnmax=c_calc(Na,2,1);
    return
end

%% if node is a reservoir:
ResID = sim_param.nodes.ReservoirID;
%concentration is the same as the initial concentration
if ismember(Na, ResID(:))
    cnmin=c_calc(Na,1,1);
    cnmax=c_calc(Na,2,1);
    return
end

%% if node is a tank:
TankID = sim_param.tanks.ID;
if ismember(Na, TankID(:))
    [cnmin, cnmax]= Algorithm_4( Na, kc, sim_param );
    return
end

%% Find all pipes that bring water into node Na:
Ain = sim_param.IncidenceMat;
Qu = sim_param.links.FlowUpper;
Ql = sim_param.links.FlowLower;
pin=[];
pin1=find(Ain(Na,:)==-1); %according to convention
for l=pin1
    if (Qu(kc,l)>0 || Ql(kc,l)>0)
        pin=[pin; l];
    end
end
pin2=find(Ain(Na,:)==1); %opposite to convention
for l=pin2
    if (Qu(kc,l)<0 || Ql(kc,l)<0)
        pin=[pin; -l];
    end
end

%% Min-Max concentration of pipes that bring water into Na:
i=1;
clmin=zeros(1,length(pin));
clmax=zeros(1,length(pin));
for l = pin'
    [clmin(i), clmax(i)] = Algorithm_2(Na,abs(l),kc,sim_param);
    i=i+1;
end

%% If only one pipe brings water into node Na:
if (size(clmin)==[1 1])
    cnmin=clmin;
    cnmax=clmax;
    c_calc(Na,1,kc)=cnmin;
    c_calc(Na,2,kc)=cnmax;
    return
end

%% if no water comes into the node:
if isempty(pin)
Klower = sim_param.tanks.DecayRateLower(1); %assume node decay same as tank decay
Kupper = sim_param.tanks.DecayRateUpper(1);
tq = sim_param.time.QualityStep;
if (~isnan(c_calc(Na,1,kc-1)) && ~isnan(c_calc(Na,2,kc-1)))
    cnmin = c_calc(Na,1,kc-1)*exp(Klower*tq);
    cnmax = c_calc(Na,2,kc-1)*exp(Kupper*tq);
    return
else
    [ c_calc(Na,1,kc-1), c_calc(Na,2,kc-1) ] = Algorithm_5(Na, kc-1, sim_param);
    cnmin = c_calc(Na,1,kc-1)*exp(Klower*tq);
    cnmax = c_calc(Na,2,kc-1)*exp(Kupper*tq);
    return
end
end

%% Pipe junction Min-Max concentration 
[cnmin, cnmax] = Algorithm_3([pin clmin' clmax'],kc, Na, sim_param);
c_calc(Na,1,kc)=cnmin;
c_calc(Na,2,kc)=cnmax;

end


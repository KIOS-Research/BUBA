function [ctmin, ctmax] = Algorithm_4( Na, kc, sim_param )
% ALGORITHM_4
% Finds chlorine concentration bounds for the water inside a storage tank

% Calculates chlorine concentration in tank at time kc
% Given that: 
%   Na is a tank node
%   Tank cannot be draining and filling at the same time  
%   Only one pipe connected to tank
%   Assumption that Ql and Qu have the same sign

%% Extract function parameters:
global c_calc
tq = sim_param.time.QualityStep;
Ain = sim_param.IncidenceMat;
Qu = sim_param.links.FlowUpper;
Ql = sim_param.links.FlowLower;
TankID = sim_param.tanks.ID;
TankV = sim_param.tanks.Volume;
KlTank = sim_param.tanks.DecayRateLower;
KuTank = sim_param.tanks.DecayRateUpper;

%% Find if tank is filling or draining:
TankIndex=find(TankID==Na);
pin=find(Ain(Na,:)==-1); %according to initial convention it is filling
init_filling=1;
if isempty(pin)
    pin=find(Ain(Na,:)==1); %according to initial convention it is draining
    init_filling=0;
end
if (init_filling && Ql(kc-1,pin)>0) || (~init_filling && Ql(kc-1,pin)<0)
    filling =1;
else
    filling =0;
end

%% Calculate previous concentration if not available
if (isnan(c_calc(Na,1,kc-1)) || isnan(c_calc(Na,2,kc-1)))
    [c_calc(Na,1,kc-1), c_calc(Na,2,kc-1)] = Algorithm_4(Na, kc-1, sim_param);
end

%% IF tank is filling:
if filling
    [cinmin, cinmax ] = Algorithm_2( Na, pin, kc-1, sim_param );
    Qinmin = min(abs(Ql(kc-1,pin)),abs(Qu(kc-1,pin)));
    Qinmax = max(abs(Ql(kc-1,pin)),abs(Qu(kc-1,pin)));
    
    ctmin = (cinmin*Qinmin)/(-KlTank(Na)*TankV(kc,TankIndex)) * (1-exp(KlTank(Na)*tq))+ ...
            (c_calc(Na,1,kc-1)*TankV(kc-1,TankIndex)/TankV(kc,TankIndex))* exp(KlTank(Na)*tq);
    ctmax = (cinmax*Qinmax)/(-KuTank(Na)*TankV(kc,TankIndex)) * (1-exp(KuTank(Na)*tq))+ ...
            (c_calc(Na,2,kc-1)*TankV(kc-1,TankIndex)/TankV(kc,TankIndex))* exp(KuTank(Na)*tq);

%% IF tank is draining:
else
   ctmin = c_calc(Na,1,kc-1)*exp(KlTank(Na)*tq);
   ctmax = c_calc(Na,2,kc-1)*exp(KuTank(Na)*tq); 
end

%%
c_calc(Na,1,kc) = ctmin;
c_calc(Na,2,kc) = ctmax;
end

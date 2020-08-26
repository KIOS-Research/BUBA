function [ cpmin, cpmax ] = Algorithm_2( Na, p, kc, sim_param )
%ALGORITHM_2 
% Finds the interval of possible chlorine contributions of a single pipe to a node

% Finds the minimum concentration cpmin that pipe p can contribute to 
% node Na at time instant kc.

%% Extract function parameters:
global c_calc
tq = sim_param.time.QualityStep;
Kl = sim_param.links.DecayRateLower;
Ku = sim_param.links.DecayRateUpper;

%% If initial conditions are reached
if kc<=1
    cpmin=c_calc(Na,1,1);
    cpmax=c_calc(Na,2,1);
    return
end

%% Find possible upstream node Naus or Na and water detention times
[ Naus, delay ] = Algorithm_1( Na, p, kc, sim_param );

%% Find concentration for every Upstream node Delay pair:
cpmin=10; % cp = infinity
cpmax=0;
i=1;
for d=delay'
%find the concentration for the upstream node at step kc-delay
        if (isnan(c_calc(Naus(i),1,kc-d)) || isnan(c_calc(Naus(i),2,kc-d)))
            [cNausmin, cNausmax] = Algorithm_5(Naus(i),kc-d,sim_param);
        else
            cNausmin = c_calc(Naus(i),1,kc-d);
            cNausmax = c_calc(Naus(i),2,kc-d);                       
        end
%calculate concentration for current delay and upstream node:
        cpminnew = exp(Kl(p)*tq*(d))*cNausmin; %%% d+1 for time uncertainty
        cpmaxnew = exp(Ku(p)*tq*(d-1))*cNausmax; %%% d-1 for time uncertainty
%save min-max concentration:
        cpmin = min(cpmin,cpminnew);
        cpmax = max(cpmax,cpmaxnew);    
        i=i+1;
end
    
end



function [ Nup, delay] = Algorithm_1( Na, p, kc, sim_param )
%ALGORITHM_1 
% Finds the upstream node Nup and detention time (delay) of a water parcel moving in a single pipe

% Finds the upstream node Nup of node Na through pipe p at each time step
% in which the flows add up
%Given that:
%  We have the same hydraulic and quality step

%% Extract function parameters:
Qu = sim_param.links.FlowUpper;
Ql = sim_param.links.FlowLower;
tq = sim_param.time.QualityStep;
L = sim_param.links.Length;
Ar = sim_param.links.Area;
Ain = sim_param.IncidenceMat;

%% Find from which node the water in pipe p flows from by initial convention
Naus=find(Ain(:,p)==1); %Node Active Up-Stream

%% If Naus is the same as the active node, then use the node at the other  
%  end of pipe and convert flows to negative
if ( Naus==Na )
    Naus=find(Ain(:,p)==-1); 
    Nupl=NaN; Nupu=NaN;
    xl=0; xu=0; d=1;
    Nup=[]; delay=[];
    
    while ( Nupl~=Nupu )
        xl=xl+tq*(-Ql(kc-d,p))/Ar(p);
        xu=xu+tq*(-Qu(kc-d,p))/Ar(p);
        
        if (xl<0); Nupl=Na; Nup=[Nup; Na]; delay=[delay; d]; 
        else if (xl>L(p)); Nupl=Naus; Nup=[Nup; Naus]; delay=[delay; d];
        else Nupl=NaN;
        end
        end
        
        if (xu<0); Nupu=Na; Nup=[Nup; Na]; delay=[delay; d];
        else if (xu>L(p)); Nupu=Naus; Nup=[Nup; Naus]; delay=[delay; d];
        else Nupu=NaN; 
        end
        end        
            
        d=d+1;
        
        %terminate if initial time step is reached
        if (kc-d)<=1
            Nup=[Nup; Naus]; Nupl=Naus; Nupu=Naus; delay=[delay; kc-1];
        end
        
    end
    
%% If Naus is not the same as Na then use positive flows    
else
    Nupl=NaN; Nupu=NaN;
    xl=0; xu=0; d=1;
    Nup=[]; delay=[];
    
    while ( Nupl~=Nupu )
        xl=xl+tq*Ql(kc-d,p)/Ar(p);
        xu=xu+tq*Qu(kc-d,p)/Ar(p);
        
        if (xl<0); Nupl=Na; Nup=[Nup; Na]; delay=[delay; d]; 
        else if (xl>L(p)); Nupl=Naus; Nup=[Nup; Naus]; delay=[delay; d];
        else Nupl=NaN;
        end
        end
        
        if (xu<0); Nupu=Na; Nup=[Nup; Na]; delay=[delay; d];
        else if (xu>L(p)); Nupu=Naus; Nup=[Nup; Naus]; delay=[delay; d];
        else Nupu=NaN; 
        end
        end        
            
        d=d+1;
        
        %terminate if initial time step is reached
        if (kc-d)<=1
            Nup=[Nup; Naus]; Nupl=Naus; Nupu=Naus; delay=[delay; (kc-1)];
        end
        
    end             
end       

end


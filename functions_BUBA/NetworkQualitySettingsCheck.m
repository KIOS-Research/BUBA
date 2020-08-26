function [] = NetworkQualitySettingsCheck(inpname)
d=epanet(inpname);
errormsg=0;
disp(sprintf('')); %empty initial line

%% check node source type (setpoint or concentration)
type = d.getNodeSourceType;
if any(any(isempty(type)))
    disp(sprintf('error: Not all nodes have a Source Type.'));
    errormsg=1;
end

%% check that hydraulic step is an integer multiple of the quality step
th = d.getTimeHydraulicStep;
tq = d.getTimeQualityStep;

if (th>=tq && mod(th,tq)==0)
else
    disp(sprintf('error: Hydraulic step not an integer multiple of the quality step.'));
    errormsg=1;
end

%% check settings
if errormsg
    error('Input network error.')
end

%%
d.unload
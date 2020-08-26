function [cnmin,cnmax] = Algorithm_3(pipeinfo,kc,Na,sim_param)
% ALGORITHM_3
% Finds the chlorine concentration bounds of a node that receives water from multiple pipes
% using optimization

% pipeinfo contains: [(vector of pipe indices), (vector of each pipe
% detention time), (minimum concentration that each pipe brings)]

%% Extract function parameters:
% global  Ql Qu infeaslog
Qu = sim_param.links.FlowUpper;
Ql = sim_param.links.FlowLower;
pipe=abs(pipeinfo(:,1));
pipenb=length(pipe);
pipe_sgn=sign(pipeinfo(:,1));
cpmin=pipeinfo(:,2);
cpmax=pipeinfo(:,3);

%% Simple formulation of linear program coefficients: f,A,b,Aeq,beq
%% A, b
Ib=[eye(pipenb); -eye(pipenb)];
Qbl=[];Qbu=[];
% convert flows to positive and find range
for i=1:pipenb
    if pipe_sgn(i)==-1
        if Qu(kc,pipe(i))<0; Qlower=-Qu(kc,pipe(i));
        else Qlower=0;
        end
        if Ql(kc,pipe(i))<0; Qupper=-Ql(kc,pipe(i));
        else Qupper=0;
        end
    else %define bounds according to convention
        if Qu(kc,pipe(i))>0; Qupper=Qu(kc,pipe(i));
        else Qupper=0;
        end
        if Ql(kc,pipe(i))>0; Qlower=Ql(kc,pipe(i));
        else Qlower=0;
        end
    end
 % if a bound is negative means that at that 
    %value the pipe does not bring water into the node      
    Qbu=[Qbu; Qupper];
    Qbl=[Qbl; Qlower];
end
Qb=[Qbu; -Qbl];
A=[Ib -Qb; zeros(1,pipenb) -1];
b=zeros(2*pipenb+1,1);

%% Aeq, beq
Aeq=[ones(1,pipenb) 0];
beq=1;

%% f' minimization
if all(cpmin == cpmin(1))
    cnmin=cpmin(1);
else
%% 
f=[];
for i=1:pipenb
    f=[f cpmin(i)];
end
f=[f 0]';

%% solving the linear program: 
% http://www.mathworks.com/help/optim/ug/linprog.html
options=optimset('Display','none');
[x,fval,exitflag] = linprog(f,A,b,Aeq,beq,[],[],[],options);
cnmin=fval;
if exitflag~=1; 
    infeaslog=[infeaslog; 0 exitflag Na kc];  cnmin=(Qbl'*cpmin/sum(Qbl)); 
%     savefile = sprintf('cmin_exitflag=%.0f_Na=%.0f_kc=%.0f', exitflag, Na, kc);
%     save(['saved_simulations/' savefile '.mat']);
end % in case of linprog infeasibility
end

%% f' minimization
if all(cpmax == cpmax(1))
    cnmax=cpmax(1);
else
%% 
f=[];
for i=1:pipenb
    f=[f cpmax(i)];
end
f=-[f 0]';
%% solving the linear program: 
% http://www.mathworks.com/help/optim/ug/linprog.html
options=optimset('Display','none');
[x,fval,exitflag] = linprog(f,A,b,Aeq,beq,[],[],[],options);
cnmax=-fval;
if exitflag~=1;
    infeaslog=[infeaslog; 1 exitflag Na kc]; cnmax=(Qbu'*cpmax/sum(Qbu)); 
%     savefile = sprintf('cmax_exitflag=%.0f_Na=%.0f_kc=%.0f', exitflag, Na, kc);
%     save(['saved_simulations/' savefile '.mat']);
end % in case of linprog infeasibility
end

end


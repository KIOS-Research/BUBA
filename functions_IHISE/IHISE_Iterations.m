function [Q, h, info] = IHISE_Iterations(d,n,restankLevels,LinkStatus,Kbnd,qbnd,Qbnd,hbnd,Pcoef, info)
%IHISE: Iterative Hydraulic Interval State Estimator

%% Initialize
Kl=Kbnd(:,1);
Ku=Kbnd(:,2);
qextl=qbnd(:,1);
qextu=qbnd(:,2);
Acoef=Pcoef.Acoef;
Bcoef=Pcoef.Bcoef;
Ccoef=Pcoef.Ccoef;
pumpInd=d.LinkPumpIndex;
np=length(Qbnd);
nu=length(hbnd);
iter=2;

%% Construct Network Matrices
[ A12, A21, hext, RTL ] = Network_Matrices( d, restankLevels,LinkStatus  );
    
%% Find initial approximate interval lines
linkInd=1:length(Qbnd);
[ slope, beta ] = Interval_Lines(n, Kl, Ku, Qbnd, linkInd, Pcoef, pumpInd, iter);

%% Define Lambda, A and B matrices 
Ll = diag(slope(:,1));
Lu = diag(slope(:,2));

Aeq = [eye(np)     zeros(np,np)  A12];
Aeq=sparse(Aeq);
Beq = -hext;

Aineq =[zeros(nu,np)  A21         zeros(nu,nu);
        zeros(nu,np) -A21         zeros(nu,nu);
        -eye(np)      Ll          zeros(np,nu);
         eye(np)     -Lu          zeros(np,nu);
        zeros(1,np)   RTL         zeros(1,nu) ;
        zeros(1,np)  -RTL         zeros(1,nu) ];
Aineq=sparse(Aineq);
Bineq =[  qextu
         -qextl
         -beta(:,1)
          beta(:,2)
          sum(qextu)
         -sum(qextl)];

zbnd(:,1) = (slope(:,1).*Qbnd(:,1))+beta(:,1);
zbnd(:,2) = (slope(:,2).*Qbnd(:,2))+beta(:,2);
LB = [ zbnd(:,1);
       Qbnd(:,1);
       hbnd(:,1)];
UB = [ zbnd(:,2);
       Qbnd(:,2);
       hbnd(:,2)];    

model.A = [Aineq;Aeq];
model.rhs = [Bineq;Beq];
model.lb = LB;
model.ub = UB;
sizeIneq = size(Aineq,1);
sizeEq   = size(Aeq,1);
model.sense(1:sizeIneq) = '<';
model.sense((sizeIneq+1):(sizeIneq + sizeEq)) = '=';
model.modelsense = 'min';
params.method = 2;
params.OutputFlag = 0;
    
%% Iterations
tolerance=1;
maxIter=20;
tic
iter=2;
criterion(1)=-inf;
criterion(2)=inf;
while ((abs(criterion(iter)-criterion(iter-1)) > tolerance) && (iter < maxIter))
    iter=iter+1;
    skipCriterion = 0;
    Q=Qbnd;
    h=hbnd;    
%% Run linear program
    for i=(np+1):(2*np)
        if any(LB>UB); error('Lower bound is larger than upper bound');end
        % Minimization Pipe (i-np)
        W=zeros(1,2*np+nu);
        W(i)=1;
        model.obj = W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xl=result.x;
            Q(i-np,1)=xl(i);
            LB(i)=xl(i);
            LB(i-np)=slope(i-np,1)*xl(i)+beta(i-np,1);
            model.lb = LB;
        else
            skipCriterion = 1;
            disp(result.status)
        end
        
        % Maximization Pipe (i-np)
        model.obj = -W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xu=result.x;
            Q(i-np,2)=xu(i);
            UB(i)=xu(i);
            UB(i-np)=slope(i-np,2)*xu(i)+beta(i-np,2);
            model.ub = UB;
        else
            skipCriterion = 1;
            disp(result.status)
        end
        
        % Recalculate pipe approximate lines
        linkInd=(i-np);
        [ slope((linkInd),:), beta((linkInd),:)] = Interval_Lines(n, Kl(linkInd), Ku(linkInd), Q(linkInd,:), linkInd, Pcoef, pumpInd, iter);
        % Redifine inequality matrices
        Ll = diag(slope(:,1));
        Lu = diag(slope(:,2));
        Aineq =[zeros(nu,np)  A21         zeros(nu,nu);
                zeros(nu,np) -A21         zeros(nu,nu);
                -eye(np)      Ll          zeros(np,nu);
                 eye(np)     -Lu          zeros(np,nu);
                zeros(1,np)   RTL         zeros(1,nu) ;
                zeros(1,np)  -RTL         zeros(1,nu) ];
        Aineq=sparse(Aineq);
        Bineq =[  qextu
                 -qextl
                 -beta(:,1)
                  beta(:,2)
                  sum(qextu)
                 -sum(qextl)];
        model.A = [Aineq;Aeq];
        model.rhs = [Bineq;Beq];
    end
    
    
    for i=(2*np+1):(2*np+nu)
        if any(LB>UB); error('Lower bound is larger than upper bound');end
        W=zeros(1,2*np+nu);
        W(i)=1;        
        model.obj = W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xl=result.x;
            h(i-2*np,1)=xl(i);
            LB(i)=xl(i);
            model.lb = LB;           
        else
            skipCriterion = 1;
            disp(result.status)
        end
        
        model.obj = -W;
        result = gurobi(model,params);
        if strcmp(result.status,'OPTIMAL')
            xu=result.x;            
            h(i-2*np,2)=xu(i);
            UB(i)=xu(i);
            model.ub = UB;             
        else
            skipCriterion = 1;
            disp(result.status)
        end          
    end
        
    %% Redefine variable bounds
    Qbnd(:,1)=Q(:,1);
    Qbnd(:,2)=Q(:,2);
    hbnd(:,1)=h(:,1);
    hbnd(:,2)=h(:,2);
    
    if skipCriterion==0
    criterion(iter) = sum([Q(:,2); h(:,2)] - [Q(:,1); h(:,1)]); %1-norm of Q
    else
    criterion(iter) = criterion(iter-1)+2*tolerance;
    end

end
info.iterPerTimeStep = [info.iterPerTimeStep iter-2];

end


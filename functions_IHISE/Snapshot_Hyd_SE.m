function [Q, h] = Snapshot_Hyd_SE(d, restankLevels, LinkStatus, K, n, Pcoef, qext )

Acoef=Pcoef.Acoef;
Bcoef=Pcoef.Bcoef;
Ccoef=Pcoef.Ccoef;
[ A12, A21, hext ] = Network_Matrices( d, restankLevels, LinkStatus  );

%% Solve Hydraulics with Todini and Pilati method
% Initialize Newton iteration 
Qinit = 1*ones(d.LinkCount,1);
hinit = d.NodeElevations';
hinit([d.NodeReservoirIndex d.NodeTankIndex])=[];
Q = Qinit; h = hinit;
np=d.LinkCount;
nu=d.NodeJunctionCount;
itnum=20; %number of iterations
% Newton iterations
for iter=1:itnum
    %A11 matrix - Hazen-Williams
    A11=diag(K.*(abs(Q).^(n-1)));
    A11Q=A11*Q;
    for i = d.LinkPumpCount
        Qp=Q(d.LinkPumpIndex(i));
        if Qp>=0
            A11Q(d.LinkPumpIndex(i))=-(Acoef(i) - Bcoef(i)*Qp^Ccoef(i));
        else
            A11Q(d.LinkPumpIndex(i))=-(Acoef(i) + Bcoef(i)*abs(Qp)^Ccoef(i));
        end
%         A11Q(d.LinkPumpIndex(i))=-(Acoef(i)*Qp^2 + Bcoef(i)*Qp + Ccoef(i));
%         A11Q(d.LinkPumpIndex(i))=-(Acoef(i)*Qp^2 + Ccoef(i));
    end    
    %Cost function Dz
    Dz = double([A11Q  + A12*h + hext;  
                 A21*Q         - qext]); 
             
    %Cost function derivative J
    N=diag(n*ones(size(A11(:,1))));
    NA11=N*A11;
    for i = 1:d.LinkPumpCount
        Qp=Q(d.LinkPumpIndex(i));
        NA11(d.LinkPumpIndex,d.LinkPumpIndex)=-(-Bcoef(i)*Ccoef(i)*abs(Qp)^(Ccoef(i)-1));
%         NA11(d.LinkPumpIndex,d.LinkPumpIndex)=-(2*Acoef(i)*Qp + Bcoef(i));
%         NA11(d.LinkPumpIndex,d.LinkPumpIndex)=-(2*Acoef(i)*Qp);
    end
    J = [NA11    A12;
          A21   zeros(d.NodeJunctionCount)];   
       
    %Update with linprog
%     options=optimset('Display', 'off');
%     [Dx,~,EXITFLAG] = linprog(ones(1,size(J,1)),[],[],J,Dz,[],[],[],options);
%     if EXITFLAG==1
%     else
%         keyboard
%     end   

    %Update with gurobi
    model.obj = ones(1,size(J,1));
    Aeq = sparse(J);
    model.A = Aeq;
    model.rhs = Dz;
    model.lb = -inf(size(Dz));
    model.ub = inf(size(Dz));
    model.sense = '=';
    params.method = 1;
    params.OutputFlag = 0;
    result = gurobi(model,params);
    if strcmp(result.status,'OPTIMAL')
        Dx=result.x;        
    else
        disp(result.status)
    end
    
    %Calculate updated state
    DQ = -Dx(1:np); 
    Dh = -Dx(np+1 : np+nu);   
    Q = Q + DQ;
    h = h + Dh;
end

% for closed valve pipes, flow bounds are exactly zero
Q(find(abs(Q)<0.00001))=0;

end

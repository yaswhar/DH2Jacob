function max = maxMotor(Jv,V,theta1range,theta2range,drange)
syms theta1 theta2 l2 d real
n=size(Jv,2);
Jv=subs(Jv,l2,300);
obj=inv(Jv)*V;
max=zeros(n,1);
Q=[theta1;theta2;d];
lowbound=[theta1range(1) theta2range(1) drange(1)];
upbound=[theta1range(2) theta2range(2) drange(2)];
options=optimset('fmincon');
options.Display='off';
ms=MultiStart('FunctionTolerance',1e-2,'UseParallel',true);
gs = GlobalSearch(ms,'XTolerance',1e-3,'StartPointsToRun','bounds');
parpool('local');
parfor i=1:n
    fun=@(x) double(subs(-abs(obj(i,:)),Q,x'));
    p=createOptimProblem('fmincon','x0',lowbound,'objective',fun,'lb',lowbound,'ub',upbound,'options',options);
    [~,tmp]=run(gs,p);
    max(i)=-tmp;
end
delete(gcp);
end
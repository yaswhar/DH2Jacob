function [Jv, Jw] = DH2Jacob(DH, Q)
len = size(DH, 1);
Hh = sym('Hh', [4 4 len+1 1]);
DH=[0 0 0 0; DH];
syms thetaz dz alphx lx real
S=[0 -alphx dz; alphx 0 -thetaz; -dz thetaz 0];
H=Rot('z',thetaz)*Trans('z',dz)*Rot('x',alphx)*Trans('x',lx);
rev=[1 0 1 0];
tmp=sym(eye(4));
a=has(DH,Q);
for i = 2:len+1
    Hh(:,:,i-1)=tmp*subs(H,[thetaz dz alphx lx],DH(i-1,:));
    tmp=Hh(:,:,i-1);
    Jw(:,i-1)=simplify(any(a(i,:)&rev)*tmp(1:3,1:3)*[0;0;1]);
    if i==len+1
        Hh(:,:,i)=tmp*subs(H,[thetaz dz alphx lx],DH(i,:));
    end
end
for i = 2:len+1
    if ~any(a(i,:)&rev)
        Jv(:,i-1)=simplify(Hh(1:3,1:3,i-1)*[0;0;1]);
    else
        Jv(:,i-1)=simplify(subs(S,[thetaz dz alphx],(Hh(1:3,1:3,i-1)*[0;0;1])')*(Hh(1:3,end,4)-Hh(1:3,end,i-1)));
    end
end
end
function R = Rot(axis,angle)
% Input Validation
validateattributes(axis,{'char'},{'nonempty','scalar'},mfilename,'Axis',1);
%validateattributes(angle,{'numeric'},{'nonempty'},mfilename,'Angle',2);

axis=double(upper(axis))-88; %Simply mapping the Axis name to index using Ascii Code
r=[cosd(angle), sind(angle);sind(angle), cosd(angle)];%2-D Rotation

%Adding the embedding rotation row and column 
tmp1=[r(1:axis,:);0,0;r(axis+1:end,:)];
R=[tmp1(:,1:axis),[0;0;0],tmp1(:,axis+1:end)];
R(axis+1,axis+1)=1;

% Changing the sign of Sinus function below the "1" in the Rotation Matrix
a=rem(axis+1,3)+1;
b=rem(a,3)+1;
R(a,b)=-R(a,b);

% Embedding another row and column
R(4,4)=1;
end
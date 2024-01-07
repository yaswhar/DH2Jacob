function T=Trans(axis,distance)
% Input Validation
validateattributes(axis,{'char'},{'nonempty','scalar'},mfilename,'Axis',1);
%validateattributes(distance,{'numeric'},{'nonempty'},mfilename,'Distance',2);

%Homogenous transformation along "axis" for the amount of "distance" 
axis=upper(axis);
if axis=='X' 
    T=[1 0 0 distance; 0 1 0 0; 0 0 1 0; 0 0 0 1]; 
elseif axis=='Y'
    T=[1 0 0 0; 0 1 0 distance; 0 0 1 0; 0 0 0 1];
elseif axis=='Z' 
    T=[1 0 0 0; 0 1 0 0; 0 0 1 distance; 0 0 0 1];
end 
end
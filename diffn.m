function [ derivative ] = diffn( field, dx,direction )
%diffn( field,order, dx,direction,boundary )
% order  1  only at moment
% direction either x = [0 1] or y = [1 0]
% boundary is  periodic


fminus = circshift(field,1.0*direction);
fplus = circshift(field,-1.0*direction);



    derivative = (fplus - fminus)/(2.0*dx);


end


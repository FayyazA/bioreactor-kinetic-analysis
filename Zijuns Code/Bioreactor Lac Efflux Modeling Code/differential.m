function [out] = differential(x)
%differential function 
%get the differential of signals: latter - former

y = x(2:end); %latter 
len = length(x); 
y(len) = 0; 
out = y - x; 
out(end) = []; 

end
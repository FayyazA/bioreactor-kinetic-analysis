function [tf rf] = fermi(x0,xdata)
A = x0(1); k = x0(2); u = x0(3); d = x0(4); 
t = xdata(:,1);
tstep = t(2) - t(1);
t0 = x0(5);
aif = xdata(:,2);
f =(t>=t0).*( A.*(1./(1+(exp(k.*(t-t0-u))))));
rf = f;
% f = (t>=t0).*(A.*((1+u)./(1+u.*exp(t*k))));
a = zeros(1,floor(x0(4)))';
b = f(1:numel(f)-floor(x0(4)));
f = vertcat(a,b);
tf = conv(f,aif);
tf = tf(1:size(aif));
tf = tf.*tstep;
end
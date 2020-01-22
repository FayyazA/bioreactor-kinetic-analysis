function [M] = SimpleBloch(x,TI)
T1 = x(1);
E2 = exp(-TI/T1);
M = 1-E2;
end
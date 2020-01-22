function [T1] = getT1(R_T1_lacin, R_T1_lacex, R_T1_pyr)
% getT1 function: get the values of T1. 

T1 = struct('pyr', 0, 'lacin', 0, 'lacex', 0); 
T1.lacin = R_T1_lacin^(-1); 
T1.lacex = R_T1_lacex^(-1); 
T1.pyr = R_T1_pyr^(-1); 

end
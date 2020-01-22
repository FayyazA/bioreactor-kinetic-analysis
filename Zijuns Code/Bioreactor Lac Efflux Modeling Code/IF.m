function [input,k,theta,mag] = IF(x, input_condition, d_set)
%IF: different input functions for simulation and fitting. 
%get the pyruvate signal input with respect to time. 

input_para = xlsread('IF.xlsx'); 
k_mat = input_para(:,1);
theta_mat = input_para(:,2);
mag_mat = input_para(:,3);
k = k_mat(d_set,:); 
theta = theta_mat(d_set,:); 
mag = mag_mat(d_set,:); 
switch (input_condition) 
    case 1
        y = gampdf(x, k, theta); %gamma variate pdf - most realistic 
end 
input = mag * y; 
end
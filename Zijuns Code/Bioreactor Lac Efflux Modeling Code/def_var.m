function [d_set, t_input, index_s, index_e, t_len, T_interval, input_condition] = def_var(data,my_number,start) 
%def_var function: put all the variables definition in this function. 
%input: datasets 
%output: parameters 

d_set = my_number; %variable to different datasets 
d_size = size(data); 
t = data(:,2); %time points 
t_index = xlsread('t0.xlsx'); %variable to different datasets. %the time point that starts the pyruvate injection. 
t0 = t(t_index(d_set,:)); 
t_input = max(t-t0,0); %translation: match the input function to the S_pyr signal. 
index_s = start; %changeable with the start time point. 
index_e = 100; %changeable with the end time point. 
t_len = d_size(1); %100 
T_interval = 3; 
input_condition = 1; %changeable 

end
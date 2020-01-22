clc;
clear;

%load data 
data = xlsread('data3.xlsx'); %variable to different filenames 
my_number = 3;
%define variables 
[d_set, t_input, index_s, index_e, t_len, T_interval, input_condition] = def_var(data,my_number); 

%get the signal of pyruvate and the signals of intracellular and extracellular
%lactate 
S_pyr = data(:,3); %pyruvate_area data 
S_lac_in = data(:,5); 
S_lac_out = data(:,6);  

%get the specific time points to fit. 
S_lac_in = S_lac_in(index_s:index_e);
S_lac_out = S_lac_out(index_s:index_e); 
S_pyr = S_pyr(index_s:index_e);

%governing equation framework 
[Sin_diff_mat, Sout_diff_mat, Spyr_diff_mat, S_pyr_mat, Sin_mat, Sout_mat] = GE(S_lac_in, S_lac_out, S_pyr, T_interval, index_e);
scale = 1./mean(abs(cat(2,Sin_mat,Sout_mat,S_pyr_mat)),1); %weights 
[input,k,theta,mag] = IF(t_input, input_condition, d_set); 

%convex optimization fitting - CVX toolbox 
cvx_begin 
    variables Kpl R_T1_lacin k_inex R_T1_lacex R_T1_pyr F; 
    minimize(norm([scale(1)*(Sin_diff_mat - Kpl * S_pyr_mat - (- R_T1_lacin - k_inex) * Sin_mat)...
             ; scale(2)*(Sout_diff_mat - k_inex * Sin_mat - (-R_T1_lacex - F) * Sout_mat) ...
             ; scale(3)*(Spyr_diff_mat - (- R_T1_pyr - Kpl - F) * S_pyr_mat - input(1:end-1))])); 
    subject to 
        32^(-1) <= R_T1_lacin <= 30^(-1);  
        32^(-1) <= R_T1_lacex <= 30^(-1); 
        51^(-1) <= R_T1_pyr <= 50^(-1); 
cvx_end

%get T1 
T1 = getT1(R_T1_lacin, R_T1_lacex, R_T1_pyr);
time = t_input(1:99);
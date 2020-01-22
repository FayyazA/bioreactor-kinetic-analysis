function [Sin_diff_mat, Sout_diff_mat, Spyr_diff_mat, S_pyr_mat, Sin_mat, Sout_mat] = GE(S_lac_in, S_lac_out, S_pyr, T_interval, index_e) 
%GE function: governing equation terms generation.

Sin_diff_mat = differential(S_lac_in)./T_interval; 
Sout_diff_mat = differential(S_lac_out)./T_interval; 
Spyr_diff_mat = differential(S_pyr)./T_interval; 
S_pyr_mat = S_pyr(1:length(S_pyr)-1); %changed all these to length(vector) - 1
Sin_mat = S_lac_in(1:length(S_lac_in)-1); 
Sout_mat = S_lac_out(1:length(S_lac_out)-1); 

end
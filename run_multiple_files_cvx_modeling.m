function result = run_multiple_files_cvx_modeling()
theFiles = ["data1.xlsx","data2.xlsx","data3.xlsx"];
starting_entries = [14,12,11];
pyr_table = zeros(9,8);
lacin_table = zeros(9,10);
lacex_table = zeros(9,18);
for k = 1 : length(theFiles)
    disp(theFiles(k))
    data = xlsread(theFiles(k));
    [d_set, t_input, index_s, index_e, t_len, T_interval, input_condition] = def_var(data,k,starting_entries(k));
    [R_T1_pyr,R_T1_lacin,R_T1_lacex,Kpl,k_inex,F,S_pyr,S_lac_in,S_lac_out,k,theta,mag] = cvx_modelingsteps(data,d_set, t_input, index_s, index_e, t_len, T_interval, input_condition)
    model_exchange_dyn_new(S_pyr,S_lac_in,S_lac_out,R_T1_pyr,Kpl,F,R_T1_lacin,k_inex,R_T1_lacex,k,theta,mag)
  drawnow; % Force display to update immediately.
end
csvwrite("Pyr_summary_startingpoint.csv",pyr_table)
csvwrite("Lacin_summary_startingpoint.csv",lacin_table)
csvwrite("Lacex_summary_startingpoint.csv",lacex_table)
end
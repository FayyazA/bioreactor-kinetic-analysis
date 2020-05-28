function result = run_multiple_files_modeling_RTSC_3compart()
theFiles = ["13C_1_RTSC_FP_CA_C1Pyr_032813"];
starting_entries = [5,6,7,8,9,10,11,12,13,14,15];
%pyr_table = zeros(10,7);
lacin_table = zeros(22,6);
lacex_table = zeros(22,18);
%A_table = zeros(10,6);
%C_table = zeros(10,4);
previous_fitted_values = [1/46.334,0.0034343,0.12039,2.3562,21.803,2.55e+08;...
                            1/49.144,0.001976,0.011,5.5086,12.483,6.03e+07;...
                            1/47.043,0.003107,0.073098,5.6994,13.694,2.90e+08;...
                            1/49.786,0.016771,0.023133,4,16,180222488.6;...
                            1/48.022,0.012459,0.2117,3.84148864,17.96032478,5.18e+08;...
                            1/49.797,0.024542,0.011,3.4328,16,1.87e+08;...
                            1/51.771,0.011373,0.070219,4.3627,12.502,1.88e+08;...
                            1/46.248,0.016329,0.1974,4.8647,13.708,4.47e+08;...
                            1/49.249,0.0043175,0.090111,4,16.60482192,6.68e+08];
lb = [min(previous_fitted_values(:,1)),min(previous_fitted_values(:,2)),min(previous_fitted_values(:,3)),min(previous_fitted_values(:,4)),min(previous_fitted_values(:,5)),min(previous_fitted_values(:,6))];
ub = [max(previous_fitted_values(:,1)),max(previous_fitted_values(:,2)),max(previous_fitted_values(:,3)),max(previous_fitted_values(:,4)),max(previous_fitted_values(:,5)),max(previous_fitted_values(:,6))];
AIC_vec = zeros(16,3);
for k = 1 : length(starting_entries)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 04082020 with lactate added",theFiles(1),starting_entries(k))
    input_vector = data_set.Pyr_area';
    [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out_recover_old(input_vector,strcat(theFiles(1),"PyrRTSC_3comp"));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    [C,S_fitdyn2] = BioRx_Pyr_Lacin_recover_old_T1Lin_19(my_array,X,strcat(theFiles(1),"LacinRTSC_3comp")); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
    lacin_table(k,:) = C;
    %C_table(k,:) = round2_output;
    my_table_2 = array2table(data_set{:,[4,8,10]});
    int2 = smooth(table2array(my_table_2(:,2)));
    my_table_2(:,2) = array2table(int2);
    int3 = smooth(table2array(my_table_2(:,3)));
    my_table_2(:,3) = array2table(int3);
    my_table_2 = rows2vars(my_table_2);
    size_table_2 = size(my_table_2)
    my_array_2 = table2array(my_table_2(:,2:size_table(2)));
    [E,Sfitdyn3] = BioRx_3mets_Lac2pks_recover_old_T1Lin_19(my_array_2,X,C,strcat(theFiles(1),"LacexRTSC_3comp")); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    AIC = aic_3(my_array_2-Sfitdyn3,length(S_fit_dyn),9);
    AIC_vec(k,:) = AIC;
    lacex_table(k,:) = E;
  drawnow; % Force display to update immediately.
end
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("Lacin_summary_RTSC_3comp.xlsx",lacin_table)
xlswrite("Lacex_summary_RTSC_3comp.xlsx",lacex_table)
xlswrite("AIC9parameters_RTSC_3comp.xlsx",AIC_vec)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
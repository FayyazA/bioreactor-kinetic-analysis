function result = run_multiple_files_modeling_T1Lin_Less_Than19()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113_with_auto","13C_1_C1Pyr_uok262_102213_with_","UOK262_100813_5DiDS1mM_LBp5","UOK262B_100813_3_DiDS1mM_LBp5","UOK262_101113_8DiDS1mM_LBp5","UOK262_042514_7DiDS_LBp5", "C1Pyr_2_012114", "UOK262_092713_1_LBp5", "UOK262_082813_p5ml_min_manu", "UOK262_siRNAcntrl_103113_1_", "UOK262_103113_5_siRNAMCT4#7", "C1Pyr_1_013114_siRNA_MCT4#7"];
starting_entries = [13,6,11,9,8,8,18,17,12,15,11,2,12,9,4,19,10,8,13,14,10,8]%[14,4,11,9,8,8,17,20,7,15,7,2,12,9,4,19];%[13,6,11,9,8,8,18,17,12,15,11,2];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
%pyr_table = zeros(10,7);
lacin_table = zeros(22,6);
lacex_table = zeros(22,17);
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
for k = 1 : length(theFiles)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 04082020 with lactate added",theFiles(k),starting_entries(k))
    input_vector = data_set.Pyr_area';
    [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out_recover_old(input_vector,strcat(theFiles(k),"PyrFinalConfirmation7"));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    [C,S_fitdyn2] = BioRx_Pyr_Lacin_recover_old_T1Lin_19(my_array,X,strcat(theFiles(k),"LacinFinalConfirmation7")); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
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
    [E,Sfitdyn3] = BioRx_3mets_Lac2pks_recover_old_T1Lin_19(my_array_2,X,C,strcat(theFiles(k),"LacexFinalConfirmation7")); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    AIC = aic_3(my_array_2-Sfitdyn3,length(S_fit_dyn),7);
    AIC_vec(k,:) = AIC;
    lacex_table(k,:) = E;
  drawnow; % Force display to update immediately.
end
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("Lacin_summary_FinalConfirmation_moredata_7.xlsx",lacin_table)
xlswrite("Lacex_summary_FinalConfirmation_moredata_7.xlsx",lacex_table)
xlswrite("AIC7parameters_FinalConfirmation_moredata.xlsx",AIC_vec)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
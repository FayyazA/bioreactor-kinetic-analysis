function table_data = run_multiple_files_modeling_synthetic_data_x0_testing(data_vec,filename)
%theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113","13C_1_C1Pyr_uok262_102213"]
%starting_entries = [14,4,11,9,8,8,17,20,7,15,7,2];%[13,6,11,9,8,8,18,17,12,15,11,2];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
%pyr_table = zeros(10,7);
%lacin_table = zeros(12,6);
%lacex_table = zeros(12,18);
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
my_vec = ["R1P","kPL","FP","R1Lin","kLP","kMCT4","R1Lex","FL","kMCT1"];
table_data = zeros(1,18);
for k = 1 : length(my_vec)
    for p = 1:10
    input_vector = data_vec(1,:);
    [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out_starting_point_increments(input_vector,strcat(filename,"PyrFinalConfirmation"),thing(p,k),my_vec(k));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    my_array = data_vec(1:2,:);
    my_array(2,:) = smooth(data_vec(2,:));
    [C,S_fitdyn2] = BioRx_Pyr_Lacin_starting_point_increments(my_array,X,strcat(filename,"LacinFinalConfirmation"),thing(p,k),my_vec(k)); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
    lacin_table(k,:) = C;
    %C_table(k,:) = round2_output;
    my_array_2 = data_vec;
    my_array_2(2,:) = smooth(data_vec(2,:));
    my_array_2(3,:) = smooth(data_vec(3,:));
    [E,Sfitdyn3] = BioRx_3mets_Lac2pks_starting_point_increments(my_array_2,X,C,strcat(filename,"LacexFinalConfirmation"),thing(p,k),my_vec(k)); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    table_data(k,:) = E;
    drawnow; % Force display to update immediately.
    end
end
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
%xlswrite("Lacin_summary_FINALsmooth_beforestartingconditions.xlsx",lacin_table)
%xlswrite("Lacex_summary_FINALsmooth_beforestartingconditions.xlsx",lacex_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
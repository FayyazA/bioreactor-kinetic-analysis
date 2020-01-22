function result = run_multiple_files_modeling_Dids()
theFiles = ["UOK262_100813_5DiDS1mM_LBp5","UOK262B_100813_3_DiDS1mM_LBp5","UOK262_101113_8DiDS1mM_LBp5","UOK262_042514_7DiDS_LBp5"];
starting_entries = [12,9,4,19];
%pyr_table = zeros(10,7);
lacin_table = zeros(10,6);
lacex_table = zeros(10,18);
%A_table = zeros(10,6);
%C_table = zeros(10,4);
%lb = [min(previous_fitted_values(:,1)),min(previous_fitted_values(:,2)),min(previous_fitted_values(:,3)),min(previous_fitted_values(:,4)),min(previous_fitted_values(:,5)),min(previous_fitted_values(:,6))];
%ub = [max(previous_fitted_values(:,1)),max(previous_fitted_values(:,2)),max(previous_fitted_values(:,3)),max(previous_fitted_values(:,4)),max(previous_fitted_values(:,5)),max(previous_fitted_values(:,6))];
for k = 1 : length(theFiles)
    data_set = importfile("C:\Users\Owner\Downloads\UOK262_DiDs_BioRx_data",theFiles(k),starting_entries(k))
    input_vector = data_set.Pyr_area';
    [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out(input_vector,strcat(theFiles(k),"PyrDids"));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    [C,S_fitdyn2] = BioRx_Pyr_Lacin(my_array,X,strcat(theFiles(k),"LacinDids")); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
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
    [E,Sfitdyn3] = BioRx_3mets_Lac2pks(my_array_2,X,C,strcat(theFiles(k),"LacexDids")); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    lacex_table(k,:) = E;
  drawnow; % Force display to update immediately.
end
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("Lacin_summaryDids.xlsx",lacin_table)
xlswrite("Lacex_summaryDids.xlsx",lacex_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
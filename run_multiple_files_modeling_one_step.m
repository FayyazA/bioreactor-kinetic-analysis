function result = run_multiple_files_modeling_one_step()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset"]
starting_entries = [13,5,11,9,8,7,17,16,16,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
lacex_table = zeros(10,13);
for k = 1 : length(theFiles)
    data_set = importfile("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 09152017",theFiles(k),starting_entries(k))
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_table_2 = array2table(data_set{:,[4,8,10]});
    int2 = smooth(table2array(my_table_2(:,2)));
    my_table_2(:,2) = array2table(int2);
    int3 = smooth(table2array(my_table_2(:,3)));
    my_table_2(:,3) = array2table(int3);
    my_table_2 = rows2vars(my_table_2);
    size_table_2 = size(my_table_2);
    my_array_2 = table2array(my_table_2(:,2:size_table(2)));
    [fit3,results3] = BioRx_kinetics_Lac2pks_one_step(my_array_2,strcat(theFiles(k),"Lacexcs"));
    lacex_table(k,:) = results3;
  drawnow; % Force display to update immediately.
end
xlswrite("Lacex_summary_one_step_sym.xlsx",lacex_table)
end
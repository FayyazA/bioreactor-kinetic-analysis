function result = run_multiple_files_modeling_loops()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113","13C_1_C1Pyr_uok262_102213"];%["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset"]
starting_entries = [13,6,11,9,8,8,18,17,12,15,11,2];%[13,5,11,9,6,11,17,16,12,16];
result_table = zeros(132,18)
for k = 1 :length(theFiles)
    numbers = [max(starting_entries(k)-5,2),max(starting_entries(k)-4,2),max(starting_entries(k)-3,2),max(starting_entries(k)-2,2),max(starting_entries(k)-1,2),starting_entries(k)-0,starting_entries(k)+1,starting_entries(k)+2,starting_entries(k)+3,starting_entries(k)+4,starting_entries(k)+5];
    for p = 1 : length(numbers)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 02152019 with lactate added",theFiles(k),numbers(p));
        
        input_vector = data_set.Pyr_area';
        [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out(input_vector,strcat(theFiles(k),"PyrFINAL",num2str(numbers(p))));
        my_table = array2table(data_set{:,[4,8]});
        int = smooth(table2array(my_table(:,2)));
        my_table(:,2) = array2table(int);
        my_table = rows2vars(my_table);
        size_table = size(my_table);
        my_array = table2array(my_table(:,2:size_table(2)));
        [C,S_fitdyn2] = BioRx_Pyr_Lacin(my_array,X,strcat(theFiles(k),"LacinFINAL",num2str(numbers(p)))); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
        my_table_2 = array2table(data_set{:,[4,8,10]});
        int2 = smooth(table2array(my_table_2(:,2)));
        my_table_2(:,2) = array2table(int2);
        int3 = smooth(table2array(my_table_2(:,3)));
        my_table_2(:,3) = array2table(int3);
        my_table_2 = rows2vars(my_table_2);
        size_table_2 = size(my_table_2);
        my_array_2 = table2array(my_table_2(:,2:size_table(2)));
        [E,Sfitdyn3] = BioRx_3mets_Lac2pks_optimized(my_array_2,X,C,strcat(theFiles(k),"LacexFINAL",num2str(numbers(p)))); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
        result_table(11*(k-1)+p,:) = E;
    end
  drawnow; % Force display to update immediately.
end
xlswrite("startingpointanalysisnoredpoints.xlsx",result_table);
end
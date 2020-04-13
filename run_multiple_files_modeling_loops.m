function result = run_multiple_files_modeling_loops()
theFiles = ["C1Pyr_2_012114", "UOK262_092713_1_LBp5", "UOK262_082813_p5ml_min_manu", "UOK262_siRNAcntrl_103113_1_", "UOK262_103113_5_siRNAMCT4#7", "C1Pyr_1_013114_siRNA_MCT4#7"];
starting_entries = [7,7,7,7,7,7];%[13,5,11,9,6,11,17,16,12,16];
result_table = zeros(132,18)
for k = 1 :length(theFiles)
    numbers = [max(starting_entries(k)-5,2),max(starting_entries(k)-4,2),max(starting_entries(k)-3,2),max(starting_entries(k)-2,2),max(starting_entries(k)-1,2),starting_entries(k)-0,starting_entries(k)+1,starting_entries(k)+2,starting_entries(k)+3,starting_entries(k)+4,starting_entries(k)+5];
    for p = 1 : length(numbers)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 04082020 with lactate added",theFiles(k),numbers(p));
        
        input_vector = data_set.Pyr_area';
        [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out(input_vector,strcat(theFiles(k),"PyrFINAL",num2str(numbers(p))));
        my_table = array2table(data_set{:,[4,8]});
        int = smooth(table2array(my_table(:,2)));
        my_table(:,2) = array2table(int);
        my_table = rows2vars(my_table);
        size_table = size(my_table);
        my_array = table2array(my_table(:,2:size_table(2)));
        [C,S_fitdyn2] = BioRx_Pyr_Lacin_recover_old_T1Lin_19(my_array,X,strcat(theFiles(k),"LacinFINAL",num2str(numbers(p)))); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
        my_table_2 = array2table(data_set{:,[4,8,10]});
        int2 = smooth(table2array(my_table_2(:,2)));
        my_table_2(:,2) = array2table(int2);
        int3 = smooth(table2array(my_table_2(:,3)));
        my_table_2(:,3) = array2table(int3);
        my_table_2 = rows2vars(my_table_2);
        size_table_2 = size(my_table_2);
        my_array_2 = table2array(my_table_2(:,2:size_table(2)));
        [E,Sfitdyn3] = BioRx_3mets_Lac2pks_recover_old_T1Lin_19(my_array_2,X,C,strcat(theFiles(k),"LacexFINAL",num2str(numbers(p)))); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
        result_table(11*(k-1)+p,:) = E;
    end
  drawnow; % Force display to update immediately.
end
xlswrite("startingpointanalysisnoredpoints_04120.xlsx",result_table);
end
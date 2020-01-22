function result = run_multiple_files_modeling_loops_DIDs()
theFiles = ["UOK262_042514_7DiDS_LBp5"];%["UOK262_100813_5DiDS1mM_LBp5","UOK262B_100813_3_DiDS1mM_LBp5","UOK262_101113_8DiDS1mM_LBp5","UOK262_042514_7DiDS_LBp5"]
starting_entries = [19,20,21,22,23,24,25,26,27];%[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21];
result_table = zeros(80,18)
for k = 1 : length(theFiles)
    for p = 1 : length(starting_entries)
    data_set = importfile("C:\Users\Owner\Downloads\UOK262_DiDs_BioRx_data",theFiles(k),starting_entries(p));
        
        input_vector = data_set.Pyr_area';
        [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out(input_vector,strcat(theFiles(k),"PyrDids",num2str(starting_entries(p))));
        my_table = array2table(data_set{:,[4,8]});
        int = smooth(table2array(my_table(:,2)));
        my_table(:,2) = array2table(int);
        my_table = rows2vars(my_table);
        size_table = size(my_table);
        my_array = table2array(my_table(:,2:size_table(2)));
        [C,S_fitdyn2] = BioRx_Pyr_Lacin(my_array,X,strcat(theFiles(k),"LacinDids",num2str(starting_entries(p)))); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
        my_table_2 = array2table(data_set{:,[4,8,10]});
        int2 = smooth(table2array(my_table_2(:,2)));
        my_table_2(:,2) = array2table(int2);
        int3 = smooth(table2array(my_table_2(:,3)));
        my_table_2(:,3) = array2table(int3);
        my_table_2 = rows2vars(my_table_2);
        size_table_2 = size(my_table_2);
        my_array_2 = table2array(my_table_2(:,2:size_table(2)));
        [E,Sfitdyn3] = BioRx_3mets_Lac2pks(my_array_2,X,C,strcat(theFiles(k),"LacexDids",num2str(starting_entries(p)))); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
        result_table(20*(k-1)+p,:) = E;
    end
  drawnow; % Force display to update immediately.
end
xlswrite("startingpointanalysisDidsjustlastcase.xlsx",result_table);
end
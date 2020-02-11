function result = run_multiple_files_Tissue_slices_lactate_sum()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113_with_auto","13C_1_C1Pyr_uok262_102213_with_","UOK262_100813_5DiDS1mM_LBp5","UOK262B_100813_3_DiDS1mM_LBp5","UOK262_101113_8DiDS1mM_LBp5","UOK262_042514_7DiDS_LBp5"];
starting_entries = [13,6,11,9,8,8,18,17,12,15,11,2,12,9,4,19]%[14,4,11,9,8,8,17,20,7,15,7,2,12,9,4,19];%[13,6,11,9,8,8,18,17,12,15,11,2];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12];
pyr_table = zeros(10,5);
lacin_table = zeros(12,10);
lacex_table = zeros(10,18);
%A_table = zeros(10,6);
%C_table = zeros(10,4);
R1L_vec = [1/10,1/12,1/14,1/16,1/18,1/20,1/24,1/26,1/28,1/30,1/32,1/34,1/36,1/38];
for j = 1 : length(R1L_vec)
for k = 1 : length(theFiles)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 02152019 with lactate added",theFiles(k),starting_entries(k))
    input_vector = data_set{:,4}';
    my_table = array2table(data_set{:,[4,11]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    [X,S_fit_dyn,results] = BioRx_kinetics_Pyrfit_starting_out_recover_old(input_vector,strcat(theFiles(k),"PyrTissuesmooth"));
    pyr_table(k,:) = results;
    %A_table(k,:) = round1_output;
    [C,S_fitdyn2] = Tissue_slices(my_array,X,strcat(theFiles(k),"LacTissue",string(j)),R1L_vec(j)); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
    lacin_table(k,:) = C;
    %C_table(k,:) = round2_output;
    %my_table_2 = array2table(data_set{:,[4,8,10]});
    %int2 = smooth(table2array(my_table_2(:,2)));
    %my_table_2(:,2) = array2table(int2);
    %int3 = smooth(table2array(my_table_2(:,3)));
    %my_table_2(:,3) = array2table(int3);
    %my_table_2 = rows2vars(my_table_2);
    %size_table_2 = size(my_table_2)
    %my_array_2 = table2array(my_table_2(:,2:size_table(2)));
    %[E,Sfitdyn3] = BioRx_3mets_Lac2pks(my_array_2,X,C,strcat(theFiles(k),"LacexFINALsmooth")); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    %lacex_table(k,:) = E;
  drawnow; % Force display to update immediately.
end
xlswrite(strcat("Pyr_AllBounds",string(j),".xlsx"),pyr_table)
xlswrite(strcat("Lac_AllBounds",string(j),".xlsx"),lacin_table)
end

%xlswrite("Lacex_summaryFINALsmooth.xlsx",lacex_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
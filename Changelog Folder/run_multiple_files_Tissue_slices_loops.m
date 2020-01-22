function result = run_multiple_files_Tissue_slices_loops()
theFiles = ["13C_1_C1Pyr_RTSC_RW_CA_032014","13C_C1Pyr_RTSC_RE_N_032114","13C_2_C1Pyr_RTSC_BJS_CA_080814","13C_1_C1Pyr_RTSC_RE_CA_032114"];%["RTSC_IP_CA_030813","13C_1_C1Pyr_RTSC_RTSC_JS_N_0114","13C_1_C1Pyr_RTSC_CA_022014","13C_1_RTSC_JB_N_C1Pyr_053013","13C_1_RTSC_MW_N_053013","13C_1_RTSC_KB_N_C1Pyr_061113","13C_2_C1Pyr_RTSC_LG_N_032114","13C_1_C1Pyr_RTSC_LG_CA_032114"];
starting_entries = [12,9,9,12];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
%pyr_table = zeros(10,7);
result_table = zeros(88,10)
for k = 1 : 4
    numbers = [max(starting_entries(k)-5,2),max(starting_entries(k)-4,2),max(starting_entries(k)-3,2),max(starting_entries(k)-2,2),max(starting_entries(k)-1,2),starting_entries(k)-0,starting_entries(k)+1,starting_entries(k)+2,starting_entries(k)+3,starting_entries(k)+4,starting_entries(k)+5];%numbers = [2,3,4,5,6,7,8,9,10,11,12];
    for p = 1 : length(numbers)
    data_set = importfile("C:\Users\Owner\Documents\All Tissue Slices Data",theFiles(k),numbers(p))
    input_vector = data_set{:,3}';
    my_table = array2table(data_set{:,[3,4]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out(input_vector,strcat(theFiles(k),"PyrTissuesmooth"));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    [C,S_fitdyn2] = Tissue_slices(my_array,X,strcat(theFiles(k),"LacTissue")); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
    result_table(11*(k-1)+p,:) = C;
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
end

%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("LacTissue_loops_new_ones.xlsx",result_table)
%xlswrite("Lacex_summaryFINALsmooth.xlsx",lacex_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
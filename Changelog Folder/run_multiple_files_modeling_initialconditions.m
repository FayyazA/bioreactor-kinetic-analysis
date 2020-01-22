function result = run_multiple_files_modeling_initialconditions()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113","13C_1_C1Pyr_uok262_102213"]
starting_entries = [14,4,11,9,8,8,17,20,7,15,7,2];%[13,6,11,9,8,8,18,17,12,15,11,2];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
%pyr_table = zeros(10,7);
lacin_table = zeros(12,6);
lacex_table = zeros(12,18);
%A_table = zeros(10,6);
%C_table = zeros(10,4);
starting_conditions = ...
  [0.020535708 0.004881736 0.360759309 0.036909456 2.22045E-14 0.215098729 ...
   0.026701587 0.122024545 0.139319766 2.333710919 22 4.170937143E+8;
   0.01963244 0.002423244 0.006386508 0.046454538 2.22045E-14 0.13715377 ...
   0.02671713 0.061069544 2.22048E-14 5.410981191 13.19435961 5.837196759E+7 ...
   ;
   0.020365537 0.003559898 0.231616351 0.041104645 2.22175E-14 0.100848761 ...
   0.026573724 0.061150967 2.22267E-14 2.854500877 21.32582161 4.914276841E+8 ...
   ;
   0.020021714 0.012781163 0.003242933 0.041125269 1.27686E-7 0.632068997 ...
   0.028011204 0.561627252 3.31388E-14 5.740538016 11.93298589 1.420987723E+8 ...
   ;
   0.020381091 0.015410376 0.32753645 0.046454455 2.32183E-14 0.709728441 ...
   0.027391566 0.599341428 2.27995E-14 3.201945154 20.29126927 6.290172892E+8 ...
   ;
   0.020273415 0.012686629 1E-8 0.039195836 2.2247E-14 0.402267803 0.027945006 ...
   0.193588909 2.23266E-14 4.800989963 12.42781118 1.499537508E+8;
   0.020737467 0.013000598 0.115086768 0.043982019 2.22048E-14 0.075666579 ...
   0.026685083 0.100342156 2.22049E-14 2.710516382 17.0351649 2.332033975E+8 ...
   ;
   0.020658385 0.056931212 0.502591471 0.034902309 2.22045E-14 0.341560178 ...
   0.027878148 0.913079178 2.22048E-14 2.291135084 19.8631418 6.128047633E+8 ...
   ;
   0.019628326 0.003841161 0.007017945 0.035336617 2.22051E-14 0.390028825 ...
   0.028011204 1 2.22051E-14 5.671277704 12.95187681 3.68333041E+8;
   0.019860478 0.007432255 0.020226529 0.046454545 2.22045E-14 0.75338315 ...
   0.02786558 0.83262326 2.22082E-14 3.619381608 10.49787879 8.320958983E+7 ...
   ;
   0.021208772 0.014063627 0.12093495 0.04110855 2.22045E-14 0.110975804 ...
   0.026788333 0.173377381 2.22066E-14 4.093292103 17.62652438 3.857189784E+7 ...
   ;
   0.020003571 0.08 1 0.036516953 2.22045E-14 0.386249833 0.027865404 0.901884112 ...
   2.22045E-14 5.340954045 16.52785926 1.298360595E+9];
for k = 1 : length(theFiles)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 02152019 with lactate added",theFiles(k),starting_entries(k))
    input_vector = data_set.Pyr_area';
    [X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out_intialconditions(input_vector,strcat(theFiles(k),"PyrFINALsmooth"),starting_conditions(k,:));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    [C,S_fitdyn2] = BioRx_Pyr_Lacin_intialconditions(my_array,X,strcat(theFiles(k),"LacinFINALsmooth"),starting_conditions(k,:)); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
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
    [E,Sfitdyn3] = BioRx_3mets_Lac2pks_initialconditions(my_array_2,X,C,strcat(theFiles(k),"LacexFINALsmooth"),starting_conditions(k,:)); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    lacex_table(k,:) = E;
  drawnow; % Force display to update immediately.
end
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("Lacin_summary_FINALsmooth_startingconditions.xlsx",lacin_table)
xlswrite("Lacex_summary_FINALsmooth_startingconditions.xlsx",lacex_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
function [lb_table,lb2_table,lb3_table,ub_table,ub2_table,ub3_table] = run_multiple_files_modeling_tight_Fl_bounds_recursive()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113","13C_1_C1Pyr_uok262_102213"]
starting_entries = [14,4,11,9,8,8,17,20,7,15,7,2];%[13,6,11,9,8,8,18,17,12,15,11,2];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
%pyr_table = zeros(10,7);
lacin_table = zeros(12,6);
lacex_table = zeros(12,18);
%A_table = zeros(10,6);
%C_table = zeros(10,4);
lb_table = zeros(12,6);
ub_table = zeros(12,6);
lb2_table = zeros(12,6);
ub2_table = zeros(12,6);
lb3_table = zeros(12,9);
ub3_table = zeros(12,9);
previous_fitted_values = [1/46.334,0.0034343,0.12039,2.3562,21.803,2.55e+08;...
                            1/49.144,0.001976,0.011,5.5086,12.483,6.03e+07;...
                            1/47.043,0.003107,0.073098,5.6994,13.694,2.90e+08;...
                            1/49.786,0.016771,0.023133,4,16,180222488.6;...
                            1/48.022,0.012459,0.2117,3.84148864,17.96032478,5.18e+08;...
                            1/49.797,0.024542,0.011,3.4328,16,1.87e+08;...
                            1/51.771,0.011373,0.070219,4.3627,12.502,1.88e+08;...
                            1/46.248,0.016329,0.1974,4.8647,13.708,4.47e+08;...
                            1/49.249,0.0043175,0.090111,4,16.60482192,6.68e+08];
%lb = [min(previous_fitted_values(:,1)),min(previous_fitted_values(:,2)),min(previous_fitted_values(:,3)),min(previous_fitted_values(:,4)),min(previous_fitted_values(:,5)),min(previous_fitted_values(:,6))];
%ub = [max(previous_fitted_values(:,1)),max(previous_fitted_values(:,2)),max(previous_fitted_values(:,3)),max(previous_fitted_values(:,4)),max(previous_fitted_values(:,5)),max(previous_fitted_values(:,6))];
for k = 1 : length(theFiles)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 02152019 with lactate added",theFiles(k),starting_entries(k))
    input_vector = data_set.Pyr_area';
    lb = [1/51, .0001, 10^-8,1.5,10,0.1*10^8];
    ub = [1/47,.08,1.0,7,24,14*10^8];
    [X,S_fit_dyn,results,lb,ub] = BioRx_kinetics_Pyrfit_starting_out_recursive(input_vector,strcat(theFiles(k),"PyrnoIFR"),lb,ub);
    lb_table(k,:) = lb;
    ub_table(k,:) = ub;
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    %disp(lb)
    %disp(ub)
    lb2 = [lb(1),lb(2),lb(3),1/30,10^-20, 0.001];
    ub2 = [ub(1),ub(2),ub(3),1/22,0.1, 1.0];
    [C,S_fitdyn2,lb2,ub2] = BioRx_Pyr_Lacin_recursive(my_array,X,strcat(theFiles(k),"LacinnoIFR"),lb2,ub2); %[round2_output,fit2,results2] = BioRx_kinetics_Lacin_NoT1Lin(my_array,round1_output,strcat(theFiles(k),"Lacincs"));
    lacin_table(k,:) = C;
    lb2_table(k,:) = lb2;
    ub2_table(k,:) = ub2;
    %C_table(k,:) = round2_output;
    my_table_2 = array2table(data_set{:,[4,8,10]});
    int2 = smooth(table2array(my_table_2(:,2)));
    my_table_2(:,2) = array2table(int2);
    int3 = smooth(table2array(my_table_2(:,3)));
    my_table_2(:,3) = array2table(int3);
    my_table_2 = rows2vars(my_table_2);
    size_table_2 = size(my_table_2)
    my_array_2 = table2array(my_table_2(:,2:size_table(2)));
    lb3 = [lb2(1),lb2(2),lb2(3),lb2(4)-.001,lb2(5),lb2(6),1/37.7,0.3, 10^-20];
    ub3 = [ub2(1),ub2(2),ub2(3),ub2(4)+.001,ub2(5),ub2(6),1/35.7,0.7, .9 ];
    [E,Sfitdyn3,lb3,ub3] = BioRx_3mets_Lac2pks_tight_Fl_bounds_recursive(my_array_2,X,C,strcat(theFiles(k),"LacexnoIFR"),lb3,ub3); %[fit3,results3] = BioRx_kinetics_Lac2pks_NoT1s(my_array_2,round1_output,round2_output,strcat(theFiles(k),"Lacexcs"));
    lacex_table(k,:) = E;
    lb3_table(k,:) = lb3;
    ub3_table(k,:) = ub3;
  drawnow; % Force display to update immediately.
end
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("Lacin_summary_noIF_recursive.xlsx",lacin_table)
xlswrite("Lacex_summary_noIF_recursive.xlsx",lacex_table)
xlswrite("lb_3noIF.xlsx",lb_table)
xlswrite("lb2_3noIF.xlsx",lb2_table)
xlswrite("lb3_3noIF.xlsx",lb3_table)
xlswrite("ub_3noIF.xlsx",ub_table)
xlswrite("ub2_3noIF.xlsx",ub2_table)
xlswrite("ub3_3noIF.xlsx",ub3_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
disp(lb)
disp(ub)
disp(lb2)
disp(ub2)
disp(lb3)
disp(ub3)
end
function result = run_multiple_files_Tissue_slices_hp_toolbox()
theFiles = ["RTSC_IP_CA_030813","13C_1_C1Pyr_RTSC_RTSC_JS_N_0114","13C_1_C1Pyr_RTSC_CA_022014","13C_1_RTSC_JB_N_C1Pyr_053013","13C_1_RTSC_MW_N_053013","13C_1_RTSC_KB_N_C1Pyr_061113","13C_2_C1Pyr_RTSC_LG_N_032114","13C_1_C1Pyr_RTSC_LG_CA_032114","13C_1_C1Pyr_RTSC_RW_CA_032014","13C_C1Pyr_RTSC_RE_N_032114","13C_2_C1Pyr_RTSC_BJS_CA_080814","13C_1_C1Pyr_RTSC_RE_CA_032114","13C_1_C1Pyr_RTSC_DD_N_011014","13C_1_C1Pyr_RTSC_HEB_RCC_070313","13C_2_C1Pyr_RTSC_HEB_RCC_070313","13C_3_C1Pyr_RTSC_HEB_RCC_070313","13C_6_C1Pyr_RTSC_HEB_RCC_070313"];%["RTSC_IP_CA_030813","13C_1_C1Pyr_RTSC_RTSC_JS_N_0114","13C_1_C1Pyr_RTSC_CA_022014","13C_1_RTSC_JB_N_C1Pyr_053013","13C_1_RTSC_MW_N_053013","13C_1_RTSC_KB_N_C1Pyr_061113","13C_2_C1Pyr_RTSC_LG_N_032114","13C_1_C1Pyr_RTSC_LG_CA_032114"];
starting_entries = [10,11,10,11,11,11,9,11,16,9,9,14,6,14,12,9,13];%did final starting point analysis with these values[13,5,11,9,6,11,17,16,12,16];%[13,5,3,9,2,8,13,11,19]; %[13,11,11,10,6,8,15,13,12]; %[13,11,11,12,6,11,15,13,12]
%pyr_table = zeros(10,7);
lacin_table = zeros(16,10);
lacex_table = zeros(16,18);
params_est = struct;
params_fixed = struct;
%A_table = zeros(10,6);
%C_table = zeros(10,4);
for k = 1 : length(theFiles)
    data_set = importfile("C:\Users\Owner\Documents\All Tissue Slices Data",theFiles(k),starting_entries(k))
    input_vector = data_set{:,3}';
    my_table = array2table(data_set{:,[3,4]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    %[X,S_fit_dyn] = BioRx_kinetics_Pyrfit_starting_out_recover_old(input_vector,strcat(theFiles(k),"PyrTissuesmooth"));
    %pyr_table(k,:) = results1;
    %A_table(k,:) = round1_output;
    [params_fit, Sfit, ufit, error_metrics] = fit_pyr_kinetics_and_input_Tissue_Slice(my_array,3,ones(2,size_table(2)-1) * 30 * (pi/180),params_fixed,params_est,[],true,theFiles(k));
    lacex_table(k,1) = params_fit.kPL;
    lacex_table(k,2) = params_fit.R1L;
    lacex_table(k,3) = params_fit.Rinj;
    lacex_table(k,4) = params_fit.Tarrival;
    lacex_table(k,5) = params_fit.Tbolus;
    lacex_table(k,6) = params_fit.FP;
    lacex_table(k,7) = params_fit.FL;
    lacex_table(k,8) = params_fit.kLP;
    lacex_table(k,9) = params_fit.R1P;
    lacex_table(k,10) = error_metrics.pyruvate.Rsq;
    lacex_table(k,11) = error_metrics.lactate.Rsq;
    Sfit = squeeze(Sfit);
    for i=1:2
        rmse(i)=gfit(my_array(i,:),Sfit(i,:),'4');     %NORMALIZED!
    end
    lacex_table(k,12) = rmse(1);
    lacex_table(k,13) = rmse(2);
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
%xlswrite("Pyr_summarycsmath.xlsx",pyr_table)
xlswrite("LacTissue_hptoolbox.xlsx",lacex_table)
%xlswrite("Lacex_summaryFINALsmooth.xlsx",lacex_table)
%xlswrite("Atablescsmath",A_table)
%xlswrite("Ctablescsmath",C_table)
end
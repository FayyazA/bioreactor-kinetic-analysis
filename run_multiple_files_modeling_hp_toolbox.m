function result = run_multiple_files_modeling_hp_toolbox()
theFiles = ["HK2_1_062714_LBp5","HK2_2_062714_LBp5","HK2B_1_062714_LBp5","UMRC6_1_C1Pyr_062714_LBp5","UMRC6_2_C1Pyr_062714_LBp5","UMRC6B_1_C1Pyr_062714_LBp5","UOK262_100813_1_LBp5","UOK262_100813B_1_LBp5","UOK262_042514_6_LBp5","NewUOKDataset","13C_3_uok262_101113_with_auto","13C_1_C1Pyr_uok262_102213_with_","UOK262_100813_5DiDS1mM_LBp5","UOK262B_100813_3_DiDS1mM_LBp5","UOK262_101113_8DiDS1mM_LBp5","UOK262_042514_7DiDS_LBp5", "C1Pyr_2_012114", "UOK262_092713_1_LBp5", "UOK262_082813_p5ml_min_manu", "UOK262_siRNAcntrl_103113_1_", "UOK262_103113_5_siRNAMCT4#7", "C1Pyr_1_013114_siRNA_MCT4#7"];
lacex_table = zeros(22,18);
AIC_vec = zeros(22,3);
params_fixed = struct;
%params_fixed.R1P = 1/50;
%params_fixed.R1Lex = 1/36.7;
params_est = struct;
for k = 1 : length(theFiles)
    data_set = importfilenew("C:\Users\Owner\Downloads\UOK262_UMRC6_HK2_BioRx_2peaks_control_data 04082020 with lactate added",theFiles(k),1);
    input_vector = data_set.Pyr_area';
    my_table = array2table(data_set{:,[4,8]});
    int = smooth(table2array(my_table(:,2)));
    my_table(:,2) = array2table(int);
    my_table = rows2vars(my_table);
    size_table = size(my_table);
    my_array = table2array(my_table(:,2:size_table(2)));
    my_table_2 = array2table(data_set{:,[4,8,10]});
    int1 = smooth(table2array(my_table_2(:,1)));
    int2 = smooth(table2array(my_table_2(:,2)));
    my_table_2(:,1) = array2table(int1);
    my_table_2(:,2) = array2table(int2);
    int3 = smooth(table2array(my_table_2(:,3)));
    my_table_2(:,3) = array2table(int3);
    my_table_2 = rows2vars(my_table_2);
    size_table_2 = size(my_table_2);
    my_array_2 = table2array(my_table_2(:,2:size_table(2)));
    [params_fit, Sfit, ufit, error_metrics] = fit_pyr_kinetics_and_input_Lacin_Lacex(my_array_2,3,ones(3,size_table(2)-1) * 30 * (pi/180),params_fixed,params_est,[],true,theFiles(k));
    Sfit = squeeze(Sfit);
    for i=1:3
        rmse(i)=gfit(my_array_2(i,:),Sfit(i,:),'4');     %NORMALIZED!
    end
    time=(1:size(Sfit,2))*3;
    figure
    subplot(3,1,1)
    plot(time,(my_array_2(1,:)),'k*--',time,(Sfit(1,:)),'r-'); 
        legend('pyr magn data','pyr magn fit');
    subplot(3,1,2)
    plot(time,(my_array_2(2,:)),'k*--',time,(Sfit(2,:)),'r-'); 
        legend('lacin magn data','lac(in) magn fit');
    subplot(3,1,3)
    plot(time,(my_array_2(3,:)),'k*--',time,(Sfit(3,:)),'r-'); 
        legend('lacex magn data','lac(ex) magn fit');
    print(gcf,'-dtiff','-r0',strcat(theFiles(k),"allthree",".tif"))
    lacex_table(k,1) = params_fit.kPL;
    lacex_table(k,2) = params_fit.kLinLex;
    lacex_table(k,3) = params_fit.R1L;
    lacex_table(k,4) = params_fit.Rinj;
    lacex_table(k,5) = params_fit.Tarrival;
    lacex_table(k,6) = params_fit.Tbolus;
    lacex_table(k,7) = params_fit.FP;
    lacex_table(k,8) = params_fit.FL;
    lacex_table(k,9) = params_fit.kLinflux;
    lacex_table(k,10) = params_fit.kLP;
    lacex_table(k,11) = error_metrics.pyruvate.Rsq;
    lacex_table(k,12) = error_metrics.lactateintra.Rsq;
    lacex_table(k,13) = error_metrics.lactateextra.Rsq;
    lacex_table(k,14) = rmse(1);
    lacex_table(k,15) = rmse(2);
    lacex_table(k,16) = rmse(3);
    lacex_table(k,17) = params_fit.R1P;
    lacex_table(k,18) = params_fit.R1Lex;
    a = reshape(Sfit,[3,size_table(2)-1]);
    AIC = aic_3(my_array_2-a,length(a),10);
    AIC_vec(k,:) = AIC;
  drawnow; % Force display to update immediately.
end
xlswrite("Lacex_summary_hptoolbox_moredata.xlsx",lacex_table)
xlswrite("AIC7parameters_hptoolbox_moredata.xlsx",AIC_vec)
end
function [lb,ub,result_vector] = changing_boundary_conditions(S_dyn,rmse_pyr,rmse_lacin,rmse_lacex,lb,ub,result_vector)
while abs(result_vector(4) - lb(4)) > lb(4) * 10^-10 | abs(result_vector(4) - ub(4)) > ub(4) * 10^-10
    for i = 1:length(result_vector)
        if abs(result_vector(i) - lb(i)) > lb(i) * 10^-10
            lb(i) = lb(i) - 0.01 * lb(i)
        end
        if abs(result_vector(i) - ub(i)) > ub(i) * 10^-10
            ub(i) = ub(i) + 0.01 * lb(i)
        end    
        [Sfit_dyn,result_vector,rmse_pyr,rmse_lacin,rmse_lacex] = BioRx_kinetics_Lac2pks_one_step(S_dyn,"Blank",lb,ub)
        [lb,ub,result_vector] = changing_boundary_conditions(S_dyn,rmse_pyr,rmse_lacin,rmse_lacex,lb,ub,result_vector)
    end
end
end
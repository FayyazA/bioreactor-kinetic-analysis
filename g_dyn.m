function res_dyn = g_dyn(x)
        res_dyn = model_exchange_dyn(x) - S_dyn; 
        res_dyn = res_dyn(:);
end
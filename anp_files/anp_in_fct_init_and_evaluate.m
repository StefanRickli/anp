function z = anp_in_fct_init_and_evaluate(t,in_params)
    if any(strcmp(who('global'),'debug'))
        global debug;
    else
        debug = false;
    end
    
    in_data.im_pz_sorted = anp_in_fct_sort_pz(in_params);
    in_data.im_pz_combinations = anp_in_fct_pz_combinations(in_params,in_data);
    
    [in_data.interval_list,...
     in_params.actual_separation_pole,...
     in_params.actual_separation_zero] = anp_in_fct_init_interval_list(in_params,in_data);
    
    in_data.interval_list = anp_in_fct_fill_interval_list_q_to_z(in_params,in_data);
    in_data.interval_list = anp_in_fct_fill_interval_list_t_to_q(in_params,in_data);
    
    z = evaluate_interval_t(in_data.interval_list,t);
    
    if debug
        figure;
        scatter(real(z),imag(z));
        axis equal;
    end
end

function z = evaluate_interval_t(interval_list,t)
    ii_interval = 1;
    z = zeros(size(t));
    for ii = 1:length(t)
        while t(ii) > interval_list(ii_interval).t(2)
            ii_interval = ii_interval + 1;
        end
        T = interval_list(ii_interval).density_fct_handle(t(ii));
        z(ii) = interval_list(ii_interval).input_fct_handle(T);
    end
end

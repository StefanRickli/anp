function [interval_list,pole_zero_combinations,separation_pole,separation_zero] = test_init_interval_list(poles,zeros,imaginary_PolesZeros_sorted,separation_max_pole,separation_max_zero)
    pole_zero_combinations = test_pole_zero_combinations(imaginary_PolesZeros_sorted,separation_max_pole,separation_max_zero);
    
    [separation_pole,separation_zero] = determine_separations(pole_zero_combinations,separation_max_pole,separation_max_zero);
    
    im_poles = [imaginary_PolesZeros_sorted([imaginary_PolesZeros_sorted.pole] == 1).value];%poles(real(poles) == 0);
    im_zeros = [imaginary_PolesZeros_sorted([imaginary_PolesZeros_sorted.zero] == 1).value];%zeros(real(zeros) == 0);
    
    pz_interval_border_on_origin = any(  [(im_poles+separation_pole)==0 , (im_poles-separation_pole)==0 , (im_zeros+separation_zero)==0 , (im_zeros+separation_zero)==0]  );
    n_intervals = 5 + 2*length([im_poles,im_zeros]) - pz_interval_border_on_origin;
    
    interval_list = repmat(struct('type',[],'t',[NaN,NaN],'t_len',NaN,'q',[NaN,NaN],'q_len',NaN,'density_fct_handle',[],'input_fct_handle',[]),n_intervals,1);
end


function [separation_pole,separation_zero] = determine_separations(pole_zero_combinations,separation_max_pole,separation_max_zero)
    separation_margin = 0.05;   % how much free space between the nearest neighboring poles/zeros? ==> avoid that the nearest pole/zero-detours could have no straight part between them. p.31
    
    s_poles = [pole_zero_combinations(:).separation_pole];
    s_zeros = [pole_zero_combinations(:).separation_zero];
    nonzero_separation_poles = s_poles(s_poles > 0);
    nonzero_separation_zeros = s_zeros(s_zeros > 0);
    
    separation_pole = min([(1-separation_margin)*nonzero_separation_poles,separation_max_pole]);
    separation_zero = min([(1-separation_margin)*nonzero_separation_zeros,separation_max_zero]);
end
function interval_list = test_init_interval_list(poles,zeros,separation_max_pole,separation_max_zero)
%     poles = [2-3i,10i,3i,-1i,0.25i];
%     zeros = [-12,2i,-5i];
%     
%     separation_max_pole = 1/4;
%     separation_max_zero = 1/20;
    
    imaginary_PolesZeros_sorted = test_struct_sort(poles,zeros);
    pole_zero_combinations = test_pole_zero_combinations(imaginary_PolesZeros_sorted,separation_max_pole,separation_max_zero);
    
    [separation_pole,separation_zero] = determine_separations(pole_zero_combinations,separation_max_pole,separation_max_zero);
    
    im_poles = poles(real(poles) == 0);
    im_zeros = zeros(real(zeros) == 0);
    im_poles = imag(im_poles);
    im_zeros = imag(im_zeros);
    
    pz_interval_border_on_origin = any(  [(im_poles+separation_pole)==0 , (im_poles-separation_pole)==0 , (im_zeros+separation_zero)==0 , (im_zeros+separation_zero)==0]  );
    n_intervals = 5 + 2*length(imaginary_PolesZeros_sorted) - pz_interval_border_on_origin;
    
    interval_list = repmat(struct('start',NaN,'end',NaN,'density_fct_handle',[],'density_fct_arguments',[],'input_fct_handle',[],'input_fct_arguments',[]),n_intervals,1);
end


function [separation_pole,separation_zero] = determine_separations(pole_zero_combinations,separation_max_pole,separation_max_zero)
    separation_margin = 0.05;   % how much free space between the nearest neighboring poles/zeros? ==> avoid that the nearest pole/zero-detours could have no straight part between them. p.31
    separation_pole = nanmin([(1-separation_margin)*[pole_zero_combinations(:).separation_pole],separation_max_pole]);
    separation_zero = nanmin([(1-separation_margin)*[pole_zero_combinations(:).separation_zero],separation_max_zero]);
end
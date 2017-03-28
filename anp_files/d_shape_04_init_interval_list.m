function [interval_list,actual_halfsecant_pole,actual_halfsecant_zero] = ...
         d_shape_04_init_interval_list(in_params,in_data)
    % Creates a list of intervals that each belong to one single shape
    % function (i.e. a straight line or a circle), based on location of
    % poles and zeros on the imaginary axis and the requested radii of the
    % main half-circle and the detours.
    % 
    % PRE:  - A sorted list of purely imaginary poles and zeros.
    %       - A list of combinations of poles and zeros that we encounter
    %         along the D-contour. It should also identify the maximal
    %         allowed halfsecant per pole and zero.
    %       - Manual limits on the maximal halfsecant of pole and zero
    %         detours.
    % 
    % POST: - A zero-initialized list of D-contour intervals.
    
    halfsecant_pole_max =    in_params.halfsecant_pole_max;
    halfsecant_zero_max =    in_params.halfsecant_zero_max;
    halfsecant_margin =      in_params.halfsecant_margin;
    im_pz_sorted =           in_data.im_pz_sorted;
    pole_zero_combinations = in_data.im_pz_combinations;
    
    [actual_halfsecant_pole,actual_halfsecant_zero] = determine_halfsecants(pole_zero_combinations,...
                                                              halfsecant_pole_max,...
                                                              halfsecant_zero_max,...
                                                              halfsecant_margin);
    
    im_poles = [im_pz_sorted([im_pz_sorted.pole] == 1).value];
    im_zeros = [im_pz_sorted([im_pz_sorted.zero] == 1).value];
    
    pz_interval_border_on_origin = any( [(im_poles+actual_halfsecant_pole)==0,...
                                         (im_poles-actual_halfsecant_pole)==0,...
                                         (im_zeros+actual_halfsecant_zero)==0,...
                                         (im_zeros+actual_halfsecant_zero)==0]);
    
    n_intervals = 5 + 2*length([im_poles,im_zeros]) - pz_interval_border_on_origin;
    
    interval_list = repmat(struct('type',[],...
                                  't',[NaN,NaN],...
                                  't_len',NaN,...
                                  'q',[NaN,NaN],...
                                  'q_len',NaN,...
                                  'density_fct_handle',[],...
                                  'input_fct_handle',[])  ,  n_intervals,1);
end


function [halfsecant_pole,halfsecant_zero] = determine_halfsecants(pole_zero_combinations,...
                                                                   halfsecant_max_pole,...
                                                                   halfsecant_max_zero,...
                                                                   halfsecant_margin)
    
    s_poles = [pole_zero_combinations(:).halfsecant_pole];
    s_zeros = [pole_zero_combinations(:).halfsecant_zero];
    
    nonzero_halfsecant_poles = s_poles(s_poles > 0);
    nonzero_halfsecant_zeros = s_zeros(s_zeros > 0);
    
    halfsecant_pole = min([(1-halfsecant_margin)*nonzero_halfsecant_poles,...
                           halfsecant_max_pole]);
    halfsecant_zero = min([(1-halfsecant_margin)*nonzero_halfsecant_zeros,...
                           halfsecant_max_zero]);
end
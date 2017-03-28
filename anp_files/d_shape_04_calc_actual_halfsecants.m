function [halfsecant_pole,halfsecant_zero] = d_shape_04_calc_actual_halfsecants(in_params,in_data)
    % Finds the upper bound for the pole and zero detour halfsecants.
    % PRE:  - A sorted list of purely imaginary poles and zeros.
    %       - A list of combinations of poles and zeros that we encounter
    %         along the D-contour. It should also identify the maximal
    %         allowed halfsecant per pole and zero detour.
    %       - Manual limits on the maximal halfsecant of pole and zero
    %         detours.
    % POST: - min(manual_max_halfsecant,max_halfsecant_of_every_combination)
    
    pole_zero_combinations = in_data.im_pz_combinations;
    halfsecant_pole_max =    in_params.halfsecant_pole_max;
    halfsecant_zero_max =    in_params.halfsecant_zero_max;
    halfsecant_margin =      in_params.halfsecant_margin;
    
    s_poles = [pole_zero_combinations(:).halfsecant_pole];
    s_zeros = [pole_zero_combinations(:).halfsecant_zero];
    
    nonzero_halfsecant_poles = s_poles(s_poles > 0);
    nonzero_halfsecant_zeros = s_zeros(s_zeros > 0);
    
    halfsecant_pole = min([(1-halfsecant_margin)*nonzero_halfsecant_poles,...
                           halfsecant_pole_max]);
    halfsecant_zero = min([(1-halfsecant_margin)*nonzero_halfsecant_zeros,...
                           halfsecant_zero_max]);
end
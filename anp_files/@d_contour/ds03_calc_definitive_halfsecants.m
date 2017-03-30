function [] = ds03_calc_definitive_halfsecants(this)
    % Finds the upper bound for the pole and zero detour halfsecants.
    % PRE:  - A sorted list of purely imaginary poles and zeros.
    %       - A list of combinations of poles and zeros that we encounter
    %         along the D-contour. It should also identify the maximal
    %         allowed halfsecant per pole and zero detour.
    %       - Manual limits on the maximal halfsecant of pole and zero
    %         detours.
    % POST: - min(manual_max_halfsecant,max_halfsecant_of_every_combination)
    
    pole_zero_combinations = this.im_pz_combinations;
    halfsecant_pole_max =    this.halfsecant_pole_max;
    halfsecant_zero_max =    this.halfsecant_zero_max;
    halfsecant_margin =      this.halfsecant_margin;
    
    s_poles = [pole_zero_combinations(:).halfsecant_pole];
    s_zeros = [pole_zero_combinations(:).halfsecant_zero];
    
    nonzero_halfsecant_poles = s_poles(s_poles > 0);
    nonzero_halfsecant_zeros = s_zeros(s_zeros > 0);
    
    % These are the values we will stick with for the rest of the D-contour
    % init!
    % 'halfsecant_margin' reduces the secants by this factor and ultimately
    % leads to the desired minimum gap between detour circles.
    this.halfsecant_pole = min([(1-halfsecant_margin)*nonzero_halfsecant_poles,...
                                halfsecant_pole_max]);
    this.halfsecant_zero = min([(1-halfsecant_margin)*nonzero_halfsecant_zeros,...
                                halfsecant_zero_max]);
end
function pole_zero_combinations = d_shape_03_pz_combinations(in_params,in_data)
    % Determines which combinations of pole-zero occur in which order and
    % calculates the distance between the neighbors.
    % These distances are later used to form an upper bound on the radii we
    % can use for the detours.
    % 
    % PRE:  A list of purely imaginary poles and zeros that is sorted in
    %       the order they appear along the D-contour.
    % POST: A list containing info about which type of detour combinations
    %       are expected along the D-contour.
    %       
    %       For example, 'pz' means that there will be first a pole detour,
    %       then a straight piece of line along the imaginary axis and then
    %       a zero detour. In order to ensure the existence of said
    %       straight piece along the im-axis, we may need to impose an
    %       upper bound on the detour radii. Hence the list, such that we
    %       can later calculate the maximal detour radius.
    
    imaginary_PolesZeros_sorted =   in_data.im_pz_sorted;
    halfsecant_pole_max =           in_params.halfsecant_pole_max;
    halfsecant_zero_max =           in_params.halfsecant_zero_max;
    
    % this ratio determines which fraction of the distance between the
    % neighboring pole/zero is used for the separation of its detour
    halfsecant_max_ratio = halfsecant_pole_max/halfsecant_zero_max;
    
    % initialize the struct array to hold the information
    pole_zero_combinations = repmat(struct('type','xx','distance',NaN,'halfsecant_pole',NaN,'halfsecant_zero',NaN), length(imaginary_PolesZeros_sorted)-1, 1 );
    
    % iterate through every pairwise pole/zero combination along the
    % imaginary axis
    for ii = 2:length(imaginary_PolesZeros_sorted)
        % determine the type of pole or zero combination
        pole_zero_combinations(ii-1).type = [imaginary_PolesZeros_sorted(ii-1).type,imaginary_PolesZeros_sorted(ii).type];
        % what's the distance between them?
        pole_zero_combinations(ii-1).distance = abs(imaginary_PolesZeros_sorted(ii-1).value - imaginary_PolesZeros_sorted(ii).value);

        % calculate the separation lengths according to the separation ratio.
        switch pole_zero_combinations(ii-1).type
            case 'pp'
                pole_zero_combinations(ii-1).halfsecant_pole = pole_zero_combinations(ii-1).distance/2;
            case 'zz'
                pole_zero_combinations(ii-1).halfsecant_zero = pole_zero_combinations(ii-1).distance/2;
            case 'pz'
                pole_zero_combinations(ii-1).halfsecant_pole = pole_zero_combinations(ii-1).distance / (1 + 1/halfsecant_max_ratio);
                pole_zero_combinations(ii-1).halfsecant_zero = pole_zero_combinations(ii-1).distance / (halfsecant_max_ratio + 1);
            case 'zp'
                pole_zero_combinations(ii-1).halfsecant_pole = pole_zero_combinations(ii-1).distance / (1 + 1/halfsecant_max_ratio);
                pole_zero_combinations(ii-1).halfsecant_zero = pole_zero_combinations(ii-1).distance / (halfsecant_max_ratio + 1);
            otherwise
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end

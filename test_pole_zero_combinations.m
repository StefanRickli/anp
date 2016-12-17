function pole_zero_combinations = test_pole_zero_combinations(imaginary_PolesZeros_sorted,separation_max_pole,separation_max_zero)
%     poles = [2-3i,10i,3i,-1i];
%     zeros = [-12,2i,-5i];
%     
%     separation_max_pole = 1/2;
%     separation_max_zero = 1/10;
%     
%     imaginary_PolesZeros_sorted = test_struct_sort(poles,zeros);
    
    % this ratio determines which fraction of the distance between the
    % neighboring pole/zero is used for the separation of its detour
    separation_max_ratio = separation_max_pole/separation_max_zero;

    % initialize the struct array to hold the information
    pole_zero_combinations = repmat(struct('type','xx','distance',NaN,'separation_pole',NaN,'separation_zero',NaN), length(imaginary_PolesZeros_sorted)-1, 1 );
    
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
                pole_zero_combinations(ii-1).separation_pole = pole_zero_combinations(ii-1).distance/2;
            case 'zz'
                pole_zero_combinations(ii-1).separation_zero = pole_zero_combinations(ii-1).distance/2;
            case 'pz'
                pole_zero_combinations(ii-1).separation_pole = pole_zero_combinations(ii-1).distance / (1 + 1/separation_max_ratio);
                pole_zero_combinations(ii-1).separation_zero = pole_zero_combinations(ii-1).distance / (separation_max_ratio + 1);
            case 'zp'
                pole_zero_combinations(ii-1).separation_pole = pole_zero_combinations(ii-1).distance / (1 + 1/separation_max_ratio);
                pole_zero_combinations(ii-1).separation_zero = pole_zero_combinations(ii-1).distance / (separation_max_ratio + 1);
            otherwise
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end

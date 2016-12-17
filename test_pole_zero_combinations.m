function pole_zero_combinations = test_pole_zero_combinations(imaginary_PolesZeros_sorted,secant_max_pole,secant_max_zero)
%     poles = [2-3i,10i,3i,-1i];
%     zeros = [-12,2i,-5i];
%     
%     secant_max_pole = 1/2;
%     secant_max_zero = 1/10;
%     
%     imaginary_PolesZeros_sorted = test_struct_sort(poles,zeros);
    
    secant_max_ratio = secant_max_pole/secant_max_zero;

    pole_zero_combinations = repmat(struct('type','xx','distance',NaN,'secant_pole',NaN,'secant_zero',NaN), length(imaginary_PolesZeros_sorted)-1, 1 );
    for ii = 2:length(imaginary_PolesZeros_sorted)
        pole_zero_combinations(ii-1).type = [imaginary_PolesZeros_sorted(ii-1).type,imaginary_PolesZeros_sorted(ii).type];
        pole_zero_combinations(ii-1).distance = abs(imaginary_PolesZeros_sorted(ii-1).value - imaginary_PolesZeros_sorted(ii).value);

        switch pole_zero_combinations(ii-1).type
            case 'pp'
                pole_zero_combinations(ii-1).secant_pole = pole_zero_combinations(ii-1).distance/2;
            case 'zz'
                pole_zero_combinations(ii-1).secant_zero = pole_zero_combinations(ii-1).distance/2;
            case 'pz'
                pole_zero_combinations(ii-1).secant_pole = pole_zero_combinations(ii-1).distance / (1 + 1/secant_max_ratio);
                pole_zero_combinations(ii-1).secant_zero = pole_zero_combinations(ii-1).distance / (secant_max_ratio + 1);
            case 'zp'
                pole_zero_combinations(ii-1).secant_pole = pole_zero_combinations(ii-1).distance / (1 + 1/secant_max_ratio);
                pole_zero_combinations(ii-1).secant_zero = pole_zero_combinations(ii-1).distance / (secant_max_ratio + 1);
            otherwise
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end 
    end
end

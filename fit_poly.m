function c = fit_poly(a,b,ya,yb)
    % calculate delta_y
    
    % init coefficients for c{1} and every even c{i}
    
    % while min-test fails and high <= 16
    
        % define g_upper(x) = g_1(x) + g_high(x)
        
        % calculate integral of g_upper(x) over [a,b]
        
        % define g_low(x)
        
        % calculate integral of g_low(x) over [a,b]
        
        % calculate factor for g_low as
        % (int(g_upper) - delta_y) / int(g_low) =: d
        
        % test if min{g(x) := g_upper(x) - g_low(x)} > 0 over [a,b]
            % if failed: check if low can be increased.
                % if low == high - 2 then set low = 2 and high = high + 2
                % and start over.
    
    % end while
    
    % set the derivative coefficient as c_prime = c{1} + c{high} - d*c{low}
    
    % integrate c_prime --> c
    
    % calculate G(a) and set the last coefficient (for the integration
    % constant) to the difference to G(a) and ya.
    
    % return the coefficient vector c
    
end
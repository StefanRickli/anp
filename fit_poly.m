function fitted_poly_coeffs = fit_poly(a,b,ya,yb,va,vb)
    dbg = true;
    assert(b > a && ya > 0 && yb > 0);
    
    % calculate delta_y
    delta_y = yb - ya;
    
    % init coefficients for c{1} and every even c{i}
    c = get_poly_coeffs(a,b,va,vb);
    
    if dbg
        f1 = figure;
        line([a b],[0 0],'linewidth',0.1);
        hold on;
        f2 = figure;
        hold on;
        line([a b],[0 0],'linewidth',0.1);
    end
    
    % while min-test fails and high <= 16
    low = 2;
    high = 8;
    derivative_negative = true;
    while derivative_negative
        % define g_upper(x) = g_1(x) + g_high(x)
        c_upper = [zeros(1,length(c{high})-length(c{1})) c{1}] + (high/low)^2*40*c{high};
        % calculate integral of g_upper(x) over [a,b]
        G_upper = diff(polyval(polyint(c_upper),[a b]));
        
        % calculate integral of g_low(x) over [a,b]
        G_low = diff(polyval(polyint(c{low}),[a b]));
        
        % calculate factor for g_low as
        % (int(g_upper) - delta_y) / int(g_low) =: d
        d = (G_upper - delta_y) / G_low;
        
        % test if min{g(x) := g_upper(x) - g_low(x)} > 0 over [a,b]
        c_g = c_upper - d*[zeros(1,length(c{high})-length(c{low})) c{low}];
        T = a:(b-a)/200:b;
        G = polyval(c_g,T);
        
        if dbg
            fprintf('fit_poly: current high=%d, low=%d\n',high,low);
            figure(f1);plot(T,G);
            
            fitted_poly_coeffs = polyint(c_g);
            delta_ya = ya - polyval(fitted_poly_coeffs,a);
            fitted_poly_coeffs(end) = delta_ya;
            figure(f2);plot(T,polyval(fitted_poly_coeffs,T));
        end
        
        if min(G(2:end-1)) <= 0
            % if failed: check if low can be increased.
                % if low == high - 2 then set low = 2 and high = high + 2
                % and start over.
            if low == high - 2
                low = 2;
                high = high + 2;
                if high == 14
                    error('degree of poly_high >= 14');
                end
            else
                low = low + 2;
            end
        else
            derivative_negative = false;
        end     
    end
    
    % integrate c_g --> c
    fitted_poly_coeffs = polyint(c_g);
    
    % calculate G(a) and set the last coefficient (for the integration
    % constant) to the difference to G(a) and ya.
    delta_ya = ya - polyval(fitted_poly_coeffs,a);
    fitted_poly_coeffs(end) = delta_ya;
    
    fprintf('fit_poly(%.1f,%.1f,%.1f,%.1f,%.1f,%.1f) returned a polynomial with degrees high=%d and low=%d\n',a,b,ya,yb,va,vb,high,low);
    % return the coefficient vector c
end

% calculates the coefficients of pure even polynomials (with exception of
% the linear one that is directly defined by a,b,ya and yb) by 
% the analytically derived formula (see p.17), namely
% g_n(x) := -(2/(a-b))^n * (x - (a+b)/2)^n + 1
function c = get_poly_coeffs(a,b,va,vb)
    dx_inv = 2/(a-b);
    m = -(a+b)/2;
    
    c = cell(16,1);
    c{1} = [(va-vb)/(a-b),(a*vb-b*va)/(a-b)];
    for ii = 2:2:16
        c{ii} = -dx_inv^ii .* n_choose_vector(ii,m) + [zeros(1,ii),1];
    end
end

% this is a vectorized version of the n_choose_k function to accept a
% vector for the k's
function coeffs = n_choose_vector(n,v)
    coeffs = zeros(1,n+1);
    for ii = 0:n
        coeffs(ii+1) = nchoosek(n,ii) * v^(ii);
    end
end
function z = exp_Lorentz(T,a,b,ya,c_lorentz,gamma_star,r1,r2,c1,c2)
    delta_x = b - a;
    
    if ~isnan(c_lorentz)
        Lorentz_0 = c_lorentz  *  (atan( (2*T-a-b) ./ (2*gamma_star) ) + ...
                                   atan(   delta_x / (2*gamma_star) ) + ...
                                   (4*gamma_star.*(a-T)) ./ (delta_x^(2)+4*gamma_star^2));
    else
        Lorentz_0 = zeros(size(T));            
    end

    G1 = c1 * (a*exp(-r1*b) + exp(-r1*a)/r1 - exp(-r1*T)/r1 - T*exp(-r1*b));
    G2 = c2 * (a*exp( r2*a) - exp( r2*a)/r2 + exp( r2*T)/r2 - T*exp( r2*a));
    G = G1 + G2;

    z = ya+G+Lorentz_0;

end

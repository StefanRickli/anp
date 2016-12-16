function fit_mixed_exp_lorentz
    home;
    close all;
    a = 0.5;
    b = 0.7;
    
    ratio = 2; % delta_y/derivative_area
    
    va = 20;
    vb = 10;
    
    derivative_area = (va + vb)/2 * (b - a);
    
    ya = 2;
    delta_y = derivative_area * ratio;
    yb = ya + delta_y;
    
    dt = (b-a)/1000;
    T = a:dt:b;
    
    c = get_poly_coeffs(a,b,va,vb);
    
    alpha = 1;
    r1 = find_r1(a,b,yb-ya,va,vb,alpha);
    r2 = alpha*r1;
    c1 = va / (exp(-r1*a) - exp(-r1*b));
    c2 = vb / (exp( r2*b) - exp( r2*a));
    g1 = c1 * (exp(-r1*T) - exp(-r1*b));
    g2 = c2 * (exp( r2*T) - exp( r2*a));
    g = g1 + g2;
    
    G1 = c1 * (a*exp(-r1*b) + exp(-r1*a)/r1 - exp(-r1*T)/r1 - T*exp(-r1*b));
    G2 = c2 * (a*exp( r2*a) - exp( r2*a)/r2 + exp( r2*T)/r2 - T*exp( r2*a));
    G = G1 + G2;

%     r1 = 50;
%     r2 = r1;
%     h1 = vb * (1 - (exp(r1*(b-T))-1) / (exp(r1*(b-a))-1) );
%     h2 = va * (1 - (exp(r2*(T-a))-1) / (exp(r2*(b-a))-1) );
%     h = h1 + h2;
    
    
%     f0 = figure;
%     plot(T,[h1;h2;h1+h2]);
    
    f1 = figure;
    plot(T,[g1;g2;g1+g2]);
    
    
    f2 = figure;
    plot(T,[G1+ya;G2+ya;G1+G2+ya]);

end

function r1 = find_r1(a,b,delta_y,va,vb,alpha)
    delta_x = b-a;
    G = @(t) va*(1/t - delta_x/(exp(t*delta_x)-1))  +  vb*(1/(alpha*t) - delta_x/(exp(alpha*t*delta_x)-1)) - delta_y;
    %G = @(t) (alpha*va+vb) / (alpha*t)   -   delta_x * (va/(exp(t*delta_x)-1) + vb/(exp(alpha*t*delta_x)-1))  -  delta_y;
    r1 = fsolve(G,a)
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
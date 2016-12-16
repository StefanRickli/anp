function fit_mixed_exp_lorentz
    home;
    close all;
    
    %% init
    
    a = 0.6
    b = 0.7
    
    ratio = 2; % delta_y/derivative_area (area under trapezoid a,b,va,vb)
    
    va = 100
    vb = 1000
    ya = 2
    
    derivative_area = (va + vb)/2 * (b - a);
    delta_y = derivative_area * ratio;
    yb = ya + delta_y
    delta_x = b-a;
    dt = (b-a)/1000;
    T = a:dt:b;
    
    tic
    %% Lorentzian function
    Agmax = 0.9/2 * (va+vb) * (b-a);
    if delta_y > Agmax
        gamma_star = find_gamma_star(a,b);
        h = 2*atan(delta_x/(2*gamma_star))  -  (4*delta_x*gamma_star) / (delta_x^2 + 4*gamma_star^2);

        delta_y_prime = delta_y - Agmax;
        c_lorentz = delta_y_prime/h;

        x0 = (a+b)/2;
        lorentz = (c_lorentz*gamma_star) ./ ((T-x0).^2 + gamma_star^2);
        lorentz_a = (4*c_lorentz*gamma_star) / (delta_x^2 + 4*gamma_star^2);
        lorentz_0 = lorentz - lorentz_a;

        Lorentz_0 = c_lorentz  *  (atan( (2*T-a-b) ./ (2*gamma_star) ) + ...
                                   atan(   delta_x / (2*gamma_star) ) + ...
                                   (4*gamma_star.*(a-T)) ./ (delta_x^(2)+4*gamma_star^2));
    else
        delta_y_prime = 0;
        lorentz_0 = zeros(size(T));
        Lorentz_0 = zeros(size(T));
    end
    
    %% exponential functions
    alpha = 1;
    r1 = find_r1(a,b,delta_y-delta_y_prime,va,vb,alpha);
    r2 = alpha*r1;
    c1 = va / (exp(-r1*a) - exp(-r1*b));
    c2 = vb / (exp( r2*b) - exp( r2*a));
    g1 = c1 * (exp(-r1*T) - exp(-r1*b));
    g2 = c2 * (exp( r2*T) - exp( r2*a));
    g = g1 + g2;
    
    G1 = c1 * (a*exp(-r1*b) + exp(-r1*a)/r1 - exp(-r1*T)/r1 - T*exp(-r1*b));
    G2 = c2 * (a*exp( r2*a) - exp( r2*a)/r2 + exp( r2*T)/r2 - T*exp( r2*a));
    G = G1 + G2;
    
    H = ya+G+Lorentz_0;
    toc
    
    
    f1 = figure('name','lorentz_0');
    f1.Position = [6,572,560,420];
    plot(T,lorentz_0);
    title('lorentz_0');
    
    f2 = figure('name','Lorentz_0');
    f2.Position = [6,57,560,420];
    plot(T,Lorentz_0);
    title('Lorentz_0');
      
    f3 = figure('name','g''s');
    f3.Position = [576,574,560,420];
    plot(T,[g1;g2;g1+g2]);
    title('g''s');
    
    f4 = figure('name','G''s');
    f4.Position = [573,57,560,420];
    plot(T,[G1+ya;G2+ya;G1+G2+ya]);
    title('G''s');
    
    f5 = figure('name','lorentz_0+g');
    f5.Position = [1151,571,560,420];
    plot(T,lorentz_0+g);
    title('lorentz_0+g');
    
    f6 = figure('name','ya+Lorentz_0+G');
    f6.Position = [1153,57,560,420];
    plot(T,ya+G+Lorentz_0);
    title('ya+Lorentz_0+G');
    
    dH = diff(H)/dt;
    fprintf('ya=  %.2f, yb=    %.2f, va=   %.2f, vb=     %.2f\n',ya,yb,va,vb);
    fprintf('H(1)=%.2f, H(end)=%.2f, dH(1)=%.2f, dH(end)=%.2f\n',H(1),H(end),dH(1),dH(end));
end

function gamma_star = find_gamma_star(a,b)
    dx = b-a;
    h = @(gamma) 2*atan(dx/(2*gamma))  -  (4*dx*gamma) / (dx^2 + 4*gamma^2);
    
    gamma_star = fsolve(h,0.1);
end

function r1 = find_r1(a,b,delta_y,va,vb,alpha)
    delta_x = b-a;
    G = @(t) va*(1/t - delta_x/(exp(t*delta_x)-1))  +  vb*(1/(alpha*t) - delta_x/(exp(alpha*t*delta_x)-1)) - delta_y;
    r1 = fsolve(G,a);
end


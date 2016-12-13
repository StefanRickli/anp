
function polyfit_der_analysis
    
    R = 10;
    x2 = 1;
    y_pole_low = 10
    y_pole_high = y_pole_low+1/10;
    x_pole_low = y_pole_low/R;
    x_pole_high = y_pole_high/R;
    v0 = R;
    v2 = R;
    xover_perc = 4/5;
    
    c0 = find_poly2_coefficients(x_pole_low,x_pole_high,y_pole_low,y_pole_high);
    f_line_unshifted = @(x) polyval(c0,x);
    
    v11 = 3;
    v12 = (R - v11*x2) / x2^4;
    
    f_poly4 = @(x) v11.*x + v12.*x.^4;
    
    x0_roots = roots([v12 0 0 v11 -y_pole_high]);
    x0 = x0_roots(imag(x0_roots) == 0 & real(x0_roots) > 0);
    
    x1 = (1-xover_perc)*x0 + xover_perc*x2;
    y1 = f_poly4(x1);
    
    c = find_poly4_coefficients(x0,x1,x2,y_pole_high,y1,R,v0,v2);
    f_polyfit4_unstretched = @(t) polyval(c,t);
    
    stretch = x0/(y_pole_high/R)
    x_pole_low_stretched = x_pole_low*stretch;
    x_pole_high_stretched = x_pole_high*stretch;
    c3 = find_poly2_coefficients(x_pole_low_stretched,x0,y_pole_low,y_pole_high);
    f_line_shiftStretch = @(t) polyval(c3,t);
    
    
    c_poly4_stretched = poly_der_fit(x0,x1,x2,y_pole_high,y1,R,v0/stretch,v2,4);
    f_polyfit4_stretched = @(t) polyval(c_poly4_stretched,t);
        
    c_poly_stretched = poly_positive_derivative_fit(x0,x1,x2,y_pole_high,y1,R,v0/stretch,v2);
    f_polyfit_auto_stretched = @(t) polyval(c_poly_stretched,t);

    inv = 1-xover_perc;
    c_begin_poly4 = find_poly4_coefficients(0,inv*x_pole_low_stretched,x_pole_low_stretched,0,y_pole_low/2.2,y_pole_low,v11,(y_pole_high-y_pole_low)/(x_pole_high_stretched-x_pole_low_stretched));
    f_polyfit4_begin = @(t) polyval(c_begin_poly4,t);

    c_begin_auto = poly_positive_derivative_fit(0,inv*x_pole_low_stretched,x_pole_low_stretched,0,y_pole_low/2.2,y_pole_low,v11,(y_pole_high-y_pole_low)/(x_pole_high_stretched-x_pole_low_stretched));
    f_polyfit_begin = @(t) polyval(c_begin_auto,t);
    
    
    
    t0 = y_pole_low/R:x2/1000:y_pole_high/R;
    t1 = 0:x2/1000:x2;
    t2 = x0:x2/1000:x2;
    t3 = y_pole_low/R*stretch:(x0-y_pole_low/R*stretch)/10:x0;
    t4 = 0:x_pole_low_stretched/1000:x_pole_low_stretched;
    
    hold off;
    p = plot(t1,f_poly4(t1)); hold on;
    p = plot(t0,f_line_unshifted(t0)); p.Color(4) = 0.2;
    p = plot(t0+(x0-y_pole_high/R),f_line_unshifted(t0)); p.Color(4) = 0.2;
    p = plot(t2,f_polyfit4_unstretched(t2)); p.Color(4) = 0.2;
    p = plot(t3,f_line_shiftStretch(t3));
    p = plot(t2,f_polyfit4_stretched(t2)); p.Color(4) = 0.2;
    p = plot(t2,f_polyfit_auto_stretched(t2));
    p = plot(t4,f_polyfit4_begin(t4)); p.Color(4) = 0.2;
    p = plot(t4,f_polyfit_begin(t4));
    
    legend('original poly4','original line', ...
           'shifted line','unstretched poly4fit','shiftStretch line','stretched poly4fit','stretched poly auto','start','Location','northwest');

    xlim([0 x2]);
    ylim([-0.05*R R]);
    line([x0 x0],[0 R],'LineWidth',0.1);
    line([0 x2],[y_pole_high y_pole_high],'LineWidth',0.1);
    line([x1 x1],[0 R],'LineWidth',0.1);
    line([0 x2],[y1 y1],'LineWidth',0.1);
%     hold off;
%     t = x0:0.01:x2;
%     plot(t,f_polyfit4_stretched(t))
%     hold on;
    
%     h = @(t) polyval([4*c_poly4_stretched(1) 3*c_poly4_stretched(2) 2*c_poly4_stretched(3) c_poly4_stretched(4)],t);
%     plot(t,h(t));
    
%     plot(t,f_polyfit5_stretched(t));
%     plot(t,f_polyfit6_stretched(t));
%     plot(t,f_polyfit7_stretched(t));
%     plot(t,f_polyfit8_stretched(t));
%     plot(t,f_polyfit_stretched(t));
    
%     legend('poly4','auto');%,'poly5','poly6','poly7','poly8');
    
end

function c = find_poly2_coefficients(x0,x1,y0,y1)
   c = [1,x0;1,x1]\[y0;y1];
   c = fliplr(c');
end

function c = find_poly3_coefficients(x0,x1,y0,y1,v0,v1)
    A = [1, x0,   x0^2,   x0^3;
         1, x1,   x1^2,   x1^3;
         0,  1, 2*x0,   3*x0^2;
         0,  1, 2*x1,   3*x1^2];
     
    b = [y0;
         y1;
         v0;
         v1];

    c = fliplr((A\b)');
end

% fits a 4th order polynomial into 3 given points [x0 y0],[x1 y1]
% and [x2 y2] with given derivatives v0 and v2 at start and end
function c = find_poly4_coefficients(x0,x1,x2,y0,y1,y2,v0,v2)
    A = [1, x0,   x0^2,   x0^3,   x0^4;
         1, x1,   x1^2,   x1^3,   x1^4;
         1, x2,   x2^2,   x2^3,   x2^4;
         0,  1, 2*x0,   3*x0^2, 4*x0^3;
         0,  1, 2*x2,   3*x2^2, 4*x2^3];
     
    b = [y0;
         y1;
         y2;
         v0;
         v2];

    c = fliplr((A\b)');
end

function c = poly_positive_derivative_fit(x0,x1,x2,y0,y1,y2,v0,v2)
    degree = 3;
    
    derivative_negative = true;
    while(derivative_negative)
        degree = degree + 1;
        if degree == 20
            error('Automatic polyfit failed. Maximum polynomial degree reached (21)');
        end
        
        [c,ns] = poly_der_fit(x0,x1,x2,y0,y1,y2,v0,v2,degree);
        dc1 = polyder(c);
        res = c;
        close all;
        for ii = 1:1
%             c = res + ii*norm(res)*ns;
            dc1 = polyder(c);
            dc2 = polyder(dc1);
            dc3 = polyder(dc2);
            dc4 = polyder(dc3);

            figure;
            x = x0:(x2-x0)/1000:x2;
            plot(x,[polyval(c,x);polyval(dc1,x)/10;polyval(dc2,x)/100;polyval(dc3,x)/1000;polyval(dc4,x)/10000]);
            line([x1 x1],[0 y2]);
            h = gca;
            h.XAxisLocation = 'origin';
            xlim([x0 x2]);
            %ylim(500*[-2 3]); title(num2str(ii));
            legend('c','dc1','dc2','dc3','dc4');
        end
        polyroots = roots(dc1);
        if any(imag(polyroots) == 0 & real(polyroots) > x0 & real(polyroots) < x2)
            
        else
            %c = poly_der_fit(x0,x1,x2,y0,y1,y2,v0,v2,degree+1);
            derivative_negative = false;
            fprintf('Polynomial autofit returned polynomial of degree %i\n',degree);
        end
    end
end

function [c,ns] = poly_der_fit(x0,x1,x2,y0,a1,y2,v0,v2,degree)
    a1=0;
    assert(degree >= 4);
    A = zeros(5,degree+1);
    
    ind0 = degree:-1:0;
    A(1,:) = ones(1,degree+1)*diag(x0.^ind0);
    A(2,:) = ones(1,degree+1)*diag(x2.^ind0);
    
    pd1 = polyder(ones(1,degree+1));
    ind1 = degree-1:-1:0;
    A(3,1:degree) = pd1 * diag(x0.^ind1);
    A(4,1:degree) = pd1 * diag(x2.^ind1);
    
    pd2 = polyder(pd1);
    ind2 = degree - 2:-1:0;
    A(5,1:degree-1) = pd2 * diag(x1.^ind2);
    
%     der = 2;
%     for ii = 6:degree+1
%         pd = polyder(pd);
%         ind = degree - der:-1:0;
%         
%         A(ii,1:degree+1-der) = pd * diag(x0.^ind);
%         
%         der = der + 1;
%     end

%     for ii = 0:degree
%         A(1,ii+1) = x0^(degree-ii);
%         A(2,ii+1) = x1^(degree-ii);
%         A(3,ii+1) = x2^(degree-ii);
%         if ii < degree
%             A(4,ii+1) = (degree-ii)*x0^(degree-ii-1);
%             A(5,ii+1) = (degree-ii)*x2^(degree-ii-1);
%         end
%     end
    
%     A(6,:) = [60*x0^2 24*x0 6 zeros(1,degree - 2)];
    
    b = [y0;
        y2;
        v0;
        v2;
        a1];
    
    c = A\b;
    ns = null(A);
end

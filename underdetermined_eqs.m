function underdetermined_eqs
    x0 = 0.5;
    x1 = 0.6;
    x2 = 1;
    y0 = 2;
    y1 = 2.5;
    y2 = 10;
    v0 = 1;
    v2 = 3;    
    
    c1 = find_poly_coefficients(x0,x1,x2,y0,y1,y2,v0,v2,6,[0 0]);
    f1 = @(t) polyval(c1,t);
    
    t = -1:0.001:2;
    hold off;
    plot(t,f1(t));
    xlim([x0 x2]);ylim([0 y2]);
    line([x0 x0],[0 y2]);
    line([x0 x2],[y0 y0]);
    line([x1 x1],[0 y2]);
    line([x0 x2],[y1 y1]);
    line([x2 x2],[0 y2]);
    hold on;
    
    c2 = find_poly_coefficients(x0,x1,x2,y0,y1,y2,v0,v2,5,0);
    f2 = @(t) polyval(c2,t);
    
    plot(t,f2(t));
    
    keyboard;
end


function c = find_poly_coefficients(x0,x1,x2,y0,y1,y2,v0,v2,degree,factors)
    a1 = zeros(1,degree+1);
    a2 = zeros(1,degree+1);
    a3 = zeros(1,degree+1);
    a4 = zeros(1,degree+1);
    a5 = zeros(1,degree+1);
    
    for ii = 0:degree
        a1(ii+1) = x0^(degree-ii);
        a2(ii+1) = x1^(degree-ii);
        a3(ii+1) = x2^(degree-ii);
        if ii < degree
            a4(ii+1) = (degree-ii)*x0^(degree-ii-1);
            a5(ii+1) = (degree-ii)*x2^(degree-ii-1);
        end
    end
    
    A = [a1;a2;a3;a4;a5];
    
%     A = [x0^5, x0^4, x0^3, x0^2, x0, 1;
%         x1^5, x1^4, x1^3, x1^2, x1, 1;
%         x2^5, x2^4, x2^3, x2^2, x2, 1;
%         5*x0^4,4*x0^3,3*x0^2,2*x0,1,0;
%         5*x2^4,4*x2^3,3*x2^2,2*x2,1,0];
%         
    b = [y0;
        y1;
        y2;
        v0;
        v2];
    
    ns = null(A);
    
    c = A\b + ns * factors';
end
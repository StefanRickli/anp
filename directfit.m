function directfit
    a = 1;
    b = 5;
    ya = 1;
    yb = 5;
    va = 12;
    vb = 0;
    
    dx = b - a;
    dy = yb-ya;
    vtot=va+vb;
    
    q1 = 2*dy/vtot
    q2 = dx
    
    c0=((-3*(a*(va+vb)-b*(va+vb)+2*dy))/((a-b)*(a^(2)-2*a*b+b^(2))));
    g1 = @(x) -c0*x.^2 + c0*(a+b)*x-c0*a*b;

    g1 = @(x) -c0*x.^2 + c0*(a+b)*x-c0*a*b;
    g2 = @(x) (va-vb)/(a-b).*x + (a*vb-b*va)/(a-b);
    g = @(x) g1(x) + g2(x);
    
    y0=((-(a^(4)*c0-4*a^(3)*b*c0+3*a^(2)*(b^(2)*c0+va+vb)-6*a*(b*va+ya)+6*b*ya))/(6*(a-b)));
    
    c3 = -c0/3;
    c2 = (((va-vb)/(2*(a-b)))+(((a+b)*c0)/(2)));
    c1 = (((-b*(va-vb))/(a-b))-a*b*c0+vb);
    c0 = y0;
    
    t = a:(b-a)/1000:b;
    f = @(x) polyval([c3 c2 c1 c0],x);
    
    plot(t,[g2(t);g(t);f(t)]);
    
    
end
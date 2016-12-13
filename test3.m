function test3
    vmin = 1;
    vpole = 3;
    tp1 = 0.65;
    p1l = 3;
    
    c0 = 0;
    c1 = vmin;
    c2 = 3*p1l/tp1^2 - (2*vmin+vpole)/tp1;
    c3 = (vmin+vpole)/tp1^2 - 2*p1l/tp1^3;
    
    t = 0:0.0001:0.7;
    f = polyval([c3 c2 c1 c0],t);
    g = polyval([3*c3 2*c2 c1],t);
    h = polyval([15 0 0 0 0],t);
    plot(t,[f;g;h])
end
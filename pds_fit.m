function pds_fit
    
    a = 2;
    b = 3;
    
    c2 = -1;
    c1 = (a+b);
    c0 = -a*b;
    
    t = 1:0.01:4;
    plot(t,polyval([c2 c1 c0],t));
    lin
    e([1 4],[0 0]);
    line([1 4],(((a+b)/2)^2-a*b)*[1 1]);
    
    keyboard;
end
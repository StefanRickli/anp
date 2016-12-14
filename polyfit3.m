
function polyfit3
    t1 = 0.1;
    t2 = 0.5;
    y1 = 1;
    y2 = 10;
    v1 = 0;
    v2 = 2;
    
    c0=((-(t1^(3)*(t2*v2-y2)+t1^(2)*t2*(t2*(v1-v2)+3*y2)-t1*t2^(2)*(t2*v1+3*y1)+t2^(3)*y1))/(t1^(3)-3*t1^(2)*t2+3*t1*t2^(2)-t2^(3)));
    c1=((t1^(3)*v2+t1^(2)*t2*(2*v1+v2)-t1*t2*(t2*(v1+2*v2)+6*(y1-y2))-t2^(3)*v1)/(t1^(3)-3*t1^(2)*t2+3*t1*t2^(2)-t2^(3)));
    c2=((-(t1^(2)*(v1+2*v2)+t1*(t2*(v1-v2)-3*(y1-y2))-t2*(t2*(2*v1+v2)+3*(y1-y2))))/(t1^(3)-3*t1^(2)*t2+3*t1*t2^(2)-t2^(3)));
    c3=((t1*(v1+v2)-t2*(v1+v2)-2*(y1-y2))/(t1^(3)-3*t1^(2)*t2+3*t1*t2^(2)-t2^(3)));
    
    t = t1:(t2-t1)/1000:t2;
    f = @(t) polyval([c3 c2 c1 c0],t);
    df = @(t) polyval([3*c3 2*c2 c1],t);
    
    plot(t,[f(t);max(f(t))*df(t)/max(df(t))]);
    f(t1)
    f(t2)
    df(t1)
    df(t2)
    
    
    
    
    
    
end
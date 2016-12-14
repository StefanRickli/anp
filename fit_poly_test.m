
a = 0.5;
b = 0.7;
ya = 2;
yb = 5;
va = 10;
vb = 65;

t = a:(b-a)/1000:b;

tic
c = fit_poly(a,b,ya,yb,va,vb);
toc

f = polyval(c,t);

plot(t,f);

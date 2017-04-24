
n = 400;

times = zeros(1,n);

for ii = 1:n
    x = 10000 * rand(1,ii);
    
    tic
    y = roots(x);
    times(ii) = toc;
end

plot(times)
function polytest
    
    a = 0.1;
    b = 0.2;
    ya = 1;
    yb = 3;
    
    dx_inv = 2/(a-b);
    m = -(a+b)/2;
    
    c = initPolyCoeffs(a,b,ya,yb,dx_inv,m);
    
    t = a:(b-a)/1000:b;
    plot(t,polyval(c{4},t)-polyval(c{2},t));
%     ylim(yb*[-1 1]);
%     line((a+b)/2*[1 1],yb*[-1 1],'linewidth',0.1);
    
end

function c = initPolyCoeffs(a,b,ya,yb,dx_inv,m)
    c = cell(16,1);
    c{1} = [(ya-yb)/(a-b),(a*yb-b*ya)/(a-b)];
    for ii = 2:2:16
        c{ii} = -dx_inv^ii .* binomialcoeffs(ii,m) + [zeros(1,ii),1];
    end    
end

function res = binomialcoeffs(degree,m)
    res = zeros(1,degree+1);
    for ii = 0:degree;
        %fprintf('%i*%f\n',nchoosek(degree,ii),m^(ii));
        res(ii+1) = nchoosek(degree,ii) * m^(ii);
    end
end
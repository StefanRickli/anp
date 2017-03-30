function z = circ_normal(q,R,phi_0,x0,y0)
    assert((real(y0) == 0 && imag(y0) ~= 0) || (real(y0) ~= 0 && imag(y0) == 0) || (y0 == 0));
    y0 = real(y0) + imag(y0);
    
    z = R*exp(1i*(phi_0 + q/R)) + x0 + 1i*y0;
end
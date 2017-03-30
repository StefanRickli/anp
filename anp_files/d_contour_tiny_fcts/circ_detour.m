function z = circ_detour(q,R,sec,y_pz)
    y_pz = real(y_pz) + imag(y_pz);
    
    z = R*(exp(1i*(q/R - asin(sec/(2*R)))) - sqrt(1-(sec/(2*R))^2)) + 1i*y_pz;
end
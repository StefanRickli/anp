function R = anp_calc_main_R(poles,zeros,min_angle_contribution_at_R,max_pole_separation)
    min_R = 2.5 * max_pole_separation;
    
    pz_all = [poles,zeros];
    
    R1 = imag(pz_all) - real(pz_all) * tan(deg2rad(min_angle_contribution_at_R));
    R2 = imag(pz_all) - real(pz_all) * tan(-deg2rad(min_angle_contribution_at_R));
    
    pz_im = pz_all(abs(real(pz_all)) < 100*eps);
    R = max([abs(R1),abs(R2),abs(pz_im)*1.5,min_R]);
end
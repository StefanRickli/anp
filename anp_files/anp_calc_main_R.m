function R = anp_calc_main_R(poles,zeros,min_angle_contribution_at_R,max_halfsecant,delay)
    pz_all = [poles,zeros];
    
    R1 = imag(pz_all) - real(pz_all) * tan(deg2rad(min_angle_contribution_at_R));
    R2 = imag(pz_all) - real(pz_all) * tan(-deg2rad(min_angle_contribution_at_R));
    
    pz_im = pz_all(abs(real(pz_all)) < 100*eps);
    if isempty(pz_im)
        pz_im = 0;
    end
    
    R = max([abs(R1),abs(R2),(abs(pz_im) + 3.5*max_halfsecant)*1.5]);
    
    if delay ~= 0
        R = 2*R + max(0,exp(-1.4*log(delay))) + max(0,100*log(delay)) + 15;
    end
    
    if R > 1e4
        warning('Calculation of D-shape radius returned %f. Are you really shure?\nComputation can take a considerable amount of time and memory!\nConsider aborting by closing the empty figure, otherwise press any key to continue and THEN close the figure.\n',R);
        waitforbuttonpress;
    end
end
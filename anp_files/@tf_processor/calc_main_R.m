function R = calc_main_R(~, poles, zeros, min_angle_contribution_at_R, max_halfsecant, delay) % ignored parameter is 'this'
    % Tries to calculate an appropriate radius for the main half circle
    % such that we can investigate what happens in the origin of the
    % Nyquist plot during that part of the D-contour.
    
    %                       ^ Im
    %                       |
    %                       |
    %        radius=omega:  o--
    %                      /|   \
    %                    /  |     \
    %                  /    |       \
    %                /      |        \
    %              /min_ang°|         |
    % ------------X-------------------|-----> Re
    %           pole        |         |
    %                       |        /
    %                       |       /
    %                       |     /
    %                       |   /
    %                       |--
    %
    % For every pole and zero, calculate the required value for the
    % radius R such that the angle contribution
    % is min_angle_contribution_at_R°

    pz_all = [poles,zeros];
    R1 = imag(pz_all) - real(pz_all) * tan(deg2rad(min_angle_contribution_at_R));
    R2 = imag(pz_all) - real(pz_all) * tan(-deg2rad(min_angle_contribution_at_R));
    
    pz_im = pz_all(abs(real(pz_all)) < 100*eps);
    if isempty(pz_im)
        pz_im = 0;
    end
    
    % Require that the angle contributions of ALL poles and zeros are at
    % least 'min_angle_contribution_at_R' and if there's a purely imaginary
    % p/z far outside, have R be 'reasonably' larger (here hardcoded)
    R = max(1.5*[abs(R1),abs(R2),(abs(pz_im) + 3.5*max_halfsecant)]);
    
    % Take delay into account by adding three additional encirclements
    if delay ~= 0
        delta_w = 3*2*pi/delay;
        R = R + delta_w;
    end
    
    if R > 1e4
        warning('Calculation of D-shape radius returned %f. Are you really shure?\nComputation can take a considerable amount of time and memory!\nConsider aborting by closing the empty figure, otherwise press any key to continue and THEN close the figure.\n',R);
        waitforbuttonpress;
    end
end
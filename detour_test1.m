function [] = test2()
    
    l = 1;
    phi = 45;
    
    get_shape(10,10);
    
    t = -l/2:0.01:l/2;
    plot(shifted_circle(t,l,phi));
    axis equal;
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    line(l*[-0.6 0.6],l/2*[1 1]);
    line(l*[-0.6 0.6],-l/2*[1 1]);
    xlim(l*[-0.6 0.6]);
    ylim(l*[-0.6 0.6]);

end

function z = shifted_circle(t,l,phi)
    phi = phi*pi/180;
    R = sqrt((l/2)^2/(1-tan(phi)^2/(4+tan(phi)^2)));
    x0 = sqrt(R^2 - (l/2)^2);
    
    z = sqrt(R^2-t.^2) - x0 + t*1i;
    
end

function shape = get_shape(R,phi)
    detour_interval_max = 1/2;
    detour_angle = 45;
    pz = [0,0.1i];
    
    t_min = 0;
    t_max = 1;
    animationParams.duration = 20;
    animationParams.FPS = 20;
    
    assert(t_min >= 0 && t_min <1 && t_max > 0 && t_max <= 1 && t_min < t_max);
    animationParams.N_frames = animationParams.duration * animationParams.FPS;
    t_step = (t_max-t_min)/(animationParams.N_frames-1);
    t = t_min:t_step:t_max;
    
    n_steps = length(t);
    
    phi = deg2rad(phi);
    
    % calculate values of the small round-off circles
    r = R*sin(phi)/(1+sin(phi));
    m2 = r + i*sqrt(R^2-2*R*r);
    m4 = r - i*sqrt(R^2-2*R*r);
    
    % need to know how long the parts are
    d1 = sqrt(R^2-2*R*r);
    d2 = r*(pi/2 + phi);
    d3 = R*(pi - 2*phi);
    d4 = d2;
    d5 = d1;
    
    % length of the round part (small arcs + big arc)
    d_arcs = d2 + d3 + d4;
    
    % calculate the t-values at which two adjacent regions touch
    % the only arbitrary choice here is t12. together with t45 it
    % determines the fraction of time we need for the y-axis part of the
    % D-curve.
    % the v-values are velocity-factors
    t12 = 2/5; 
    t45 = 1-t12;
    v_ave = d_arcs / (t45 - t12);
    
    t23 = t12 + d2/v_ave;
    t34 = 1-t23;
    
    v11 = 3; % again the only arbitrary choice amongst the velocity-factors. determines the speed of the animation near the origin.
    v12 = (d1-v11*t12)/t12^4;
    v51 = v11;
    v52 = v12;
        
    ts = t_min:t_step:t_max;    % '1st row': t-values
    shape = zeros(n_steps,1);   % '2nd row': shape-values
    
    region1_idx = find(ts >=t_min          & ts <= min(t_max,t12));
    region2_idx = find(ts > max(t_min,t12) & ts <= min(t_max,t23));
    region3_idx = find(ts > max(t_min,t23) & ts <= min(t_max,t34));
    region4_idx = find(ts > max(t_min,t34) & ts <= min(t_max,t45));
    region5_idx = find(ts > max(t_min,t45) & ts <=t_max);
    
    shape(region1_idx) = (v11*ts(region1_idx) + v12*ts(region1_idx).^4)*1i; % part 1
    shape(region2_idx) = m2 + r*exp(1i*(pi - (ts(region2_idx)-t12).*(pi/2 + phi)/(t23-t12))); % part 2
    shape(region3_idx) = R*exp(1i*(pi/2 - phi - (ts(region3_idx)-t23).*(pi-2*phi)/(t34-t23))); % part 3
    shape(region4_idx) = m4 + r*exp(1i*(phi - pi/2 - (ts(region4_idx)-t34).*(pi/2 + phi)/(t45-t34))); % part 4
    shape(region5_idx) = -(v51*(1-ts(region5_idx)) + v52*(1-ts(region5_idx)).^4)*1i; % part 5
    
    % get poles and zeros which lie on the imaginary axis
    imag_pz_idx = find(real(pz) == 0);
    imaginary_pz = pz(imag_pz_idx);
    
    % find minimum distance between any pair of only-imaginary poles and
    % zeros
    mask = triu(ones(length(imaginary_pz)),1);
    [pz1 pz2] = meshgrid(imaginary_pz,imaginary_pz);
    pairwise_distances = abs(pz1.*mask - pz2.*mask);
    d_min = min(nonzeros(pairwise_distances));
    
    % determine the radii for the detours 
    % 1)    based on the calculated minimum distances between imaginary
    %       poles and zeros
    % 2)    based on the provided maximal secant length
    detour_angle = detour_angle*pi/180;
    max_pairwise_detour_radius = sqrt((d_min/2)^2/(1-tan(detour_angle)^2/(4+tan(detour_angle)^2)));
    max_parametrized_detour_radius =  sqrt((detour_interval_max/2)^2/(1-tan(detour_angle)^2/(4+tan(detour_angle)^2)));
    
    % set the parameters for all of the detour circles
    secant_detour = min(d_min,detour_interval_max);
    radius_detour = min(max_pairwise_detour_radius,max_parametrized_detour_radius);
    x_shift_detour = sqrt(radius_detour^2 - (secant_detour/2)^2);
    
    for ii = 1:length(imaginary_pz)
        detour_interval_low = imag(imaginary_pz(ii))-secant_detour/2;
        detour_interval_high = imag(imaginary_pz(ii))+secant_detour/2;
        detour_point_idx = find(real(shape) == 0 & imag(shape) > detour_interval_low & imag(shape) < detour_interval_high);
        shape(detour_point_idx) = shifted_circle(shape(detour_point_idx),imaginary_pz(ii));
    end
    
    pause;
    
    function z = shifted_circle(t,y_shift_detour)
        z = sqrt(radius_detour^2-imag(t-y_shift_detour).^2) - x_shift_detour + imag(t).*1i;
    end
    
end
